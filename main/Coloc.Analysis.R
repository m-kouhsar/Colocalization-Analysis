
get.lead.snps <- function(loci.data){
  result <- vector(mode="character" , length=0)
  for (i in 1:nrow(loci.data)) {
    if((as.numeric(loci.data[i,]$END) - as.numeric(loci.data[i,]$START)) < 2){
      result <- append(result,loci.data[i,]$SNP)# Lead SNP
    }
  }
  return(data.frame(id=result))
}
##################################################################
args<-commandArgs(TRUE)

QTL.file <- args[1]
GWAS.file <- args[2]
loci.file <- args[3]
ref.genome.prefix <- args[4]
distance_ <- as.numeric(args[5])          ## distance in bp to define regions based on LD (it will be ignored if length locus > 2)
type_ <- args[6]                         ## "quant" for quantitative traits and "cc" for binary traits
ld.threshold <- as.numeric(args[7])      ## LD threshold for defining regions (it should be between 0 and 1)
out.pref <- args[8]
#################################################################
print("Input arguments:")
print(paste0("     QTL file= ",QTL.file))
print(paste0("     GWAS Summary statistic file= ",GWAS.file))
print(paste0("     Loci file= ",loci.file))
print(paste0("     Coloc type= ",type_))
print(paste0("     Distance for selecting region based on LD (if the region length < 2)= ",distance_))
print(paste0("     Reference genome binary files name (--bfile option in plink)= ",ref.genome.prefix))
print(paste0("     Output prefix= ",out.pref))
cat('\n')

print("Loading libraries...")
suppressMessages(library(coloc))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(data.table))
options(datatable.fread.datatable=FALSE)

print("Reading QTL file...")
QTLs <- vroom(QTL.file, col_types = c(.default = "c"),num_threads=16)
QTLs <- as.data.frame(QTLs)
names(QTLs) <- toupper(names(QTLs))
print("Removing dupplicated SNPs from QTLs...")
QTLs <- QTLs[!duplicated(QTLs$SNPID),] 
QTLs$MAF <- as.numeric(QTLs$MAF)
QTLs$SE <- as.numeric(QTLs$SE)
QTLs$BETA <- as.numeric(QTLs$BETA)
QTLs$N <- as.numeric(QTLs$N)
QTLs <- QTLs[(!is.na(QTLs$BETA))&(!is.na(QTLs$SE))&(!is.na(QTLs$SNPID))&(!is.na(QTLs$MAF)),]
QTLs <- QTLs[QTLs$MAF>0,]
QTLs <- QTLs[QTLs$MAF<1,]
QTLs_coloc = list(beta = as.numeric(QTLs$BETA), # Beta/expression values for allele1 of each CPG/gene
                  varbeta = QTLs$SE^2, # Variance of each beta value
                  type = type_, # quant for quantitative or cc for binary outcome
                  snp = QTLs$SNPID, # IDs of QTLs. MUST be the same as in the GWAS dataset
                  MAF = QTLs$MAF, # Frequency of allele1 of each SNP
                  N = QTLs$N) # Sample size of the associated study.In general will be the max of available sample size values
check_ <- check_dataset(QTLs_coloc) # That functions returns NULL if everything is OK for the coloc analysis
if(is.null(check_)){
  print("QTL dataset is OK")
}

print("Reading GWAS summary statistics...")
GWAS <- vroom(file=GWAS.file,col_types = c(.default = "c"),num_threads=16)
GWAS <- as.data.frame(GWAS)
names(GWAS) <- toupper(names(GWAS))
print("Removing dupplicated SNPs from GWAS summary statistics...")
GWAS <- GWAS[!duplicated(GWAS$SNP),]
GWAS$MAF <- as.numeric(GWAS$MAF)
GWAS$SE <- as.numeric(GWAS$SE)
GWAS$BETA <- as.numeric(GWAS$BETA)
GWAS$N <- as.numeric(GWAS$N)
GWAS <- GWAS[(!is.na(GWAS$BETA))&(!is.na(GWAS$SE))&(!is.na(GWAS$SNP))&(!is.na(GWAS$MAF)),]
GWAS <- GWAS[GWAS$MAF>0,]
GWAS <- GWAS[GWAS$MAF<1,]
GWAS_coloc = list(beta = GWAS$BETA, 
                  varbeta = GWAS$SE^2,
                  type = type_,
                  snp = GWAS$SNP, 
                  MAF = GWAS$MAF,
                  N = GWAS$N)
check_ <- check_dataset(GWAS_coloc)
if(is.null(check_)){
  print("GWAS dataset is OK")
}
total.common.snps <- length(QTLs$SNPID[QTLs$SNPID %in% GWAS$SNP])
if(total.common.snps > 0)
{
  print(paste("Total number of shared SNPs between QTLs and GWAS:",total.common.snps))
  print("Reading loci file...")
  loci <- vroom(file = loci.file,col_types = c(.default = "c"),num_threads=16)
  loci <- as.data.frame(loci)
  names(loci) <- toupper(names(loci))
  loci$START <- as.numeric(loci$START)
  loci$END <- as.numeric(loci$END)
  if(!("ID" %in% names(loci))){
    if("RSID" %in% names(loci)){
      loci$ID <- loci$RSID
    }else{
      loci$ID <- loci$SNP
    }
  }
  loci$START <- as.numeric(loci$START)
  loci$END <- as.numeric(loci$END)
  
  snp.list <- get.lead.snps(loci)
  
  ld <- NULL
  if(nrow(snp.list) > 0){
    if(file.exists(paste0(loci.file,".ld"))){
      print("Reading LD file...")
      ld <- fread(paste0(loci.file,".ld"),stringsAsFactors = F , header = T)
    }else{
      print("Calculating LD...")
      write.table(snp.list,file = paste0(loci.file,"lead.snp.list"),row.names = F,col.names = F,quote = F)
      plink_cmd <- paste("plink --r2 --bfile",ref.genome.prefix , "--ld-snp-list", paste0(loci.file,"lead.snp.list") ,"--ld-window-r2",ld.threshold,"--ld-window-kb",distance_,"--out",loci.file)
      exit_code <- system(plink_cmd, ignore.stdout=T,wait = T)
      if(!(exit_code>0)){
        rm <- file.remove(paste0(loci.file,".log"))
        rm <- file.remove(paste0(loci.file,"lead.snp.list"))
        ld <- fread(paste0(loci.file,".ld"),stringsAsFactors = F , header = T)
      }else{
        stop(paste("Could not calculate LD using plin"),"\n","plink command:","\n",plink_cmd,"\n","exit code:",exit_code)
      }
    } 
  }
  
  log_ <- as.data.frame(matrix(data = "",nrow = nrow(loci),ncol = 23))
  names(log_) <- c("Coloc.Type","LD.Threshold","Distance","UniqGene.QTL","UniqSNPs.QTL","UniqSNPs.GWAS","CommonSNPs.QTL.GWAS",
                   "locus.Id","locus.Start","locus.End","locus.length","UniqSNPs.locus","CommonSNPs.QTL.locus","PP.H0.abf","PP.H1.abf",
                   "PP.H2.abf","PP.H3.abf","PP.H4.abf","sum.H3.H4","SNPs","nSNPs","Genes","nGenes")
  
  result <- vector(mode = "list",length = nrow(loci))
  names(result) <- paste0("locus.",loci$ID)
  for (i in 1:nrow(loci)) {
    print(paste("Processing locus",i,":",loci$ID[i]))
    log_$nSNPs[i] = 0
    log_$PP.H0.abf[i] = 0
    log_$PP.H1.abf[i] = 0
    log_$PP.H2.abf[i] = 0
    log_$PP.H3.abf[i] = 0
    log_$PP.H4.abf[i] = 0
    log_$sum.H3.H4[i] = 0
    log_$locus.Id[i] = loci$ID[i]
    if(!is.null(ld)){
      if(loci$SNP[i] %in% ld$SNP_A){
        ld1 <- ld[ld$SNP_A==loci$SNP[i],]
        loci$START[i] <- min(ld1$BP_B)
        loci$END[i] <- max(ld1$BP_B)
      }
    }
    log_$locus.Start[i] <- loci$START[i]
    log_$locus.End[i] <- loci$END[i]
    log_$locus.length[i] <- loci$END[i] - loci$START[i]
    index <- (GWAS$CHR == loci$CHR[i]) & (as.numeric(GWAS$POS) > loci$START[i]) & (as.numeric(GWAS$POS) < loci$END[i])
    GWAS.loci <- GWAS[index , ]
    c <- length(QTLs$SNPID[QTLs$SNPID %in% GWAS.loci$SNP])
    log_$CommonSNPs.QTL.locus[i] <- c
    log_$UniqSNPs.locus[i] <- length(unique(GWAS.loci$SNP))
    if(c > 0){
      GWAS_coloc = list(beta = GWAS.loci$BETA,
                        varbeta = GWAS.loci$SE^2,
                        type = type_,
                        snp = GWAS.loci$SNP, 
                        MAF = GWAS.loci$MAF,
                        N = GWAS.loci$N)
      
      result.coloc = coloc.abf(GWAS_coloc, QTLs_coloc)
      result.QTL = QTLs[QTLs$SNPID %in% GWAS.loci$SNP,]
      result[[i]] <- list(result.coloc = result.coloc,result.QTL = result.QTL)
      temp = as.data.frame(result.coloc$summary)
      log_$nSNPs[i] = temp["nsnps",1]
      log_$PP.H0.abf[i] = temp["PP.H0.abf",1]
      log_$PP.H1.abf[i] = temp["PP.H1.abf",1]
      log_$PP.H2.abf[i] = temp["PP.H2.abf",1]
      log_$PP.H3.abf[i] = temp["PP.H3.abf",1]
      log_$PP.H4.abf[i] = temp["PP.H4.abf",1]
      log_$sum.H3.H4[i] = temp["PP.H3.abf",1] + temp["PP.H4.abf",1]
      Genes = paste(unique(result.QTL$GENE),collapse=',')
      log_$Genes[i] = Genes
      log_$nGenes[i] = length(unique(result.QTL$GENE))
      log_$SNPs[i] <- paste(unique(result.coloc$results$snp),collapse=',')
      log_$nSNPs[i] <- length(unique(result.coloc$results$snp))
      
    }
  }
  
  log_ <- log_[order(log_$sum.H3.H4 , decreasing = T),]
  log_$Coloc.Type[1] <- type_
  log_$LD.Threshold[1] <- ld.threshold
  log_$Distance[1] <- distance_
  log_$UniqGene.QTL[1] <- length(unique(QTLs$GENE))
  log_$UniqSNPs.QTL[1] <- length(QTLs$SNPID)
  log_$UniqSNPs.GWAS[1] <- length(GWAS$SNP)
  log_$CommonSNPs.QTL.GWAS[1] <- length(QTLs$SNPID[QTLs$SNPID %in% GWAS$SNP])
  
  write.csv(log_,file=paste0(out.pref,".dist.",distance_,".coloc.",type_,".csv"),row.names=F)
  save(result,file = paste0(out.pref,".dist.",distance_,".coloc.",type_,".rdat"))
}else{
  print("There is no shared SNPs between QTLs and GWAS!")
}
