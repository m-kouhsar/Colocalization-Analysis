args<-commandArgs(TRUE)

QTL.file <- args[1]
GWAS.file <- args[2]
loci.file <- args[3]
type_ <- args[4]                         ## "quant" for quantitative traits and "cc" for binary traits
distance_ <- as.numeric(args[5])          ## distance in bp to define regions based on LD (it will be ignored if length locus > 2)
use.ld <- ifelse(tolower(args[6])=="yes",T,F)
ld.threshold <- as.numeric(args[7])      ## LD threshold for defining regions (it should be between 0 and 1)
ref.genome.prefix <- args[8]
out.pref <- args[9]
sink(file = paste0(out.pref , ".log.txt"))
#################################################################
cat("Input arguments:\n")
cat("     QTL directory= ",dirname(QTL.file),"\n")
cat("     GWAS directory= ",dirname(GWAS.file),"\n")
cat("     Loci directory= ",dirname(loci.file),"\n")
cat("\n")
cat("     QTL file= ",basename(QTL.file),"\n")
cat("     GWAS Summary statistic file= ",basename(GWAS.file),"\n")
cat("     Loci file= ",basename(loci.file),"\n")
cat("\n")
cat("     Coloc type= ",type_,"\n")
cat("     Distance for selecting region based on LD (if the region length < 2)= ",distance_, " KB\n")
cat("     Using LD to define regions for lead SNPs? ",ifelse(use.ld,"Yes","No"),"\n")
cat("     LD threshold= ",ld.threshold,"\n")
cat('\n')
cat("     Reference genome binary files name (--bfile option in plink)= ",ref.genome.prefix,"\n")
cat("     Output prefix= ",out.pref,"\n")
cat('\n')
cat("#############################################################################################\n")
cat("\n")
cat("Loading libraries...\n")
suppressMessages(library(coloc))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(data.table))
options(datatable.fread.datatable=FALSE)

cat("Reading QTL file...\n")
QTLs <- vroom(QTL.file, col_types = c(.default = "c"),num_threads=16)
QTLs <- as.data.frame(QTLs)
names(QTLs) <- toupper(names(QTLs))

QTLs.col <- c("SNP","MAF","SE","BETA","N")
if(!all(QTLs.col %in% colnames(QTLs))){
   missed.col <- QTLs.col[which(!(QTLs.col %in% colnames(QTLs)))]
   stop("The following columns are missed in QTL data: ",paste(missed.col , collapse = ","))
}

cat("Removing dupplicated SNPs from QTLs...\n")
QTLs <- QTLs[!duplicated(QTLs$SNP),] 
QTLs$MAF <- as.numeric(QTLs$MAF)
QTLs$SE <- as.numeric(QTLs$SE)
QTLs$BETA <- as.numeric(QTLs$BETA)
QTLs$N <- as.numeric(QTLs$N)
QTLs <- QTLs[(!is.na(QTLs$BETA))&(!is.na(QTLs$SE))&(!is.na(QTLs$SNP))&(!is.na(QTLs$MAF)),]
QTLs <- QTLs[QTLs$MAF>0,]
QTLs <- QTLs[QTLs$MAF<1,]
QTLs_coloc = list(beta = as.numeric(QTLs$BETA), # Beta/expression values for allele1 of each CPG/gene
                  varbeta = QTLs$SE^2, # Variance of each beta value
                  type = type_, # quant for quantitative or cc for binary outcome
                  snp = QTLs$SNP, # IDs of QTLs. MUST be the same as in the GWAS dataset
                  MAF = QTLs$MAF, # Frequency of allele1 of each SNP
                  N = QTLs$N) # Sample size of the associated study.In general will be the max of available sample size values
check_ <- check_dataset(QTLs_coloc) # That functions returns NULL if everything is OK for the coloc analysis
if(is.null(check_)){
  cat("QTL dataset is OK\n")
}

cat("Reading GWAS summary statistics...\n")
GWAS <- vroom(file=GWAS.file,col_types = c(.default = "c"),num_threads=16)
GWAS <- as.data.frame(GWAS)
names(GWAS) <- toupper(names(GWAS))

GWAS.col <- c("SNP","MAF","SE","BETA","N")
if(!all(GWAS.col %in% colnames(GWAS))){
  missed.col <- GWAS.col[which(!(GWAS.col %in% colnames(GWAS)))]
  stop(paste("The following columns are missed in GWAS data: ",paste(missed.col , collapse = ",")))
}

cat("Removing dupplicated SNPs from GWAS summary statistics...\n")
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
  cat("GWAS dataset is OK\n")
}
total.common.snps <- length(intersect(QTLs$SNP , GWAS$SNP))
if(total.common.snps > 0){
  cat("Total number of shared SNPs between QTLs and GWAS: ",total.common.snps,"\n")
  cat("Reading loci file...\n")
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
  loci$length <-  loci$END - loci$START
  if(any(loci$length < 2)){
    index <- which(loci$length < 2)
    if(use.ld){
      if(!("SNP" %in% colnames(loci))){
        if(!("RSID" %in% colnames(loci))){
          stop("SNP or RSID column not fond in loci file!")
        } else {
          loci$SNP <- loci$RSID
        }
      }
      ld.file.prefix <- paste0(loci.file,".dist.",distance_,".ld.",ld.threshold)
      if(file.exists(paste0(ld.file.prefix,".ld"))){
        cat("Reading LD file...\n")
        ld <- fread(paste0(ld.file.prefix,".ld"),stringsAsFactors = F , header = T)
      }else{
        cat("Calculating LD...\n")
        write.table(loci$SNP[index],file =paste0(ld.file.prefix,".SNPList.txt"),col.names = F , row.names = F,quote = F )
        plink_cmd <- paste("plink --r2 --bfile",ref.genome.prefix , "--ld-snp-list", paste0(ld.file.prefix,".SNPList.txt") ,"--ld-window-r2",ld.threshold,"--ld-window-kb",distance_,"--out",ld.file.prefix)
        exit_code <- system(plink_cmd, ignore.stdout=T,wait = T)
        if(!(exit_code>0)){
          rm <- file.remove(paste0(ld.file.prefix,".log"))
          ld <- fread(paste0(ld.file.prefix,".ld"),stringsAsFactors = F , header = T)
        }else{
          rm <- file.remove(paste0(ld.file.prefix,".log"))
          stop("Could not calculate LD using plin!","\n","plink command:","\n",plink_cmd,"\n","exit code:",exit_code)
        }
      }
      for (i in 1:length(index)) {
        if(loci$SNP[i] %in% ld$SNP_A){
          ld1 <- ld[ld$SNP_A==loci$SNP[i],]
          loci$START[i] <- min(ld1$BP_B)
          loci$END[i] <- max(ld1$BP_B)
        }
      }
    }else{
      loci$START[index] <- loci$START[index] - (distance_*1000)
      loci$END[index] <- loci$END[index] + (distance_*1000)
    }
  }
  
  log_ <- as.data.frame(matrix(data = "",nrow = nrow(loci),ncol = 23))
  names(log_) <- c("Coloc.Type","LD.Threshold","Distance","UniqGene.QTL","UniqSNPs.QTL","UniqSNPs.GWAS","CommonSNPs.QTL.GWAS",
                   "locus.Id","locus.Start","locus.End","locus.length","UniqSNPs.locus","CommonSNPs.QTL.locus","PP.H0.abf","PP.H1.abf",
                   "PP.H2.abf","PP.H3.abf","PP.H4.abf","sum.H3.H4","SNPs","nSNPs","Genes","nGenes")
  
  result <- vector(mode = "list",length = nrow(loci))
  names(result) <- paste0("locus.",loci$ID)
  flag <- F
  cat('\n')
  for (i in 1:nrow(loci)) {
    cat("Processing locus ",i,": ",loci$ID[i],"\n")
    log_$nSNPs[i] = 0
    log_$PP.H0.abf[i] = 0
    log_$PP.H1.abf[i] = 0
    log_$PP.H2.abf[i] = 0
    log_$PP.H3.abf[i] = 0
    log_$PP.H4.abf[i] = 0
    log_$sum.H3.H4[i] = 0
    log_$locus.Id[i] = loci$ID[i]
    log_$locus.Start[i] <- loci$START[i]
    log_$locus.End[i] <- loci$END[i]
    log_$locus.length[i] <- loci$END[i] - loci$START[i]
    index <- (GWAS$CHR == loci$CHR[i]) & (as.numeric(GWAS$POS) > loci$START[i]) & (as.numeric(GWAS$POS) < loci$END[i])
    GWAS.loci <- GWAS[index , ]
	log_$UniqSNPs.locus[i] <- length(unique(GWAS.loci$SNP))
	cat("Region ",loci$ID[i]," length: ",log_$locus.length[i],"\n")
	cat("Number of unique SNPs in ",loci$ID[i]," region: ",log_$UniqSNPs.locus[i],"\n")
    c <- length(intersect(QTLs$SNP , GWAS.loci$SNP))
    cat("Number of shared SNPs between QTLs and ",loci$ID[i]," region: ",c,"\n")
    log_$CommonSNPs.QTL.locus[i] <- c
    if(c > 0){
      flag <- T
      GWAS_coloc = list(beta = GWAS.loci$BETA,
                        varbeta = GWAS.loci$SE^2,
                        type = type_,
                        snp = GWAS.loci$SNP, 
                        MAF = GWAS.loci$MAF,
                        N = GWAS.loci$N)
      
      result.coloc = coloc.abf(GWAS_coloc, QTLs_coloc)
      result.QTL = QTLs[QTLs$SNP %in% GWAS.loci$SNP,]
      result[[i]] <- list(result.coloc = result.coloc,result.QTL = result.QTL)
      temp = as.data.frame(result.coloc$summary)
      log_$nSNPs[i] = temp["nsnps",1]
      log_$PP.H0.abf[i] = temp["PP.H0.abf",1]
      log_$PP.H1.abf[i] = temp["PP.H1.abf",1]
      log_$PP.H2.abf[i] = temp["PP.H2.abf",1]
      log_$PP.H3.abf[i] = temp["PP.H3.abf",1]
      log_$PP.H4.abf[i] = temp["PP.H4.abf",1]
      log_$sum.H3.H4[i] = temp["PP.H3.abf",1] + temp["PP.H4.abf",1]
      cat("Sum of H3 and H4 probability related to ",loci$ID[i]," region: ",log_$sum.H3.H4[i],"\n")
      Genes = paste(unique(result.QTL$GENE),collapse=',')
      log_$Genes[i] = Genes
      log_$nGenes[i] = length(unique(result.QTL$GENE))
      log_$SNPs[i] <- paste(unique(result.coloc$results$snp),collapse=',')
      log_$nSNPs[i] <- length(unique(result.coloc$results$snp))
      
    }
    cat('\n')
  }
  
  log_ <- log_[order(log_$sum.H3.H4 , decreasing = T),]
  log_$Coloc.Type[1] <- type_
  log_$LD.Threshold[1] <- ld.threshold
  log_$Distance[1] <- distance_
  log_$UniqGene.QTL[1] <- length(unique(QTLs$GENE))
  log_$UniqSNPs.QTL[1] <- length(QTLs$SNP)
  log_$UniqSNPs.GWAS[1] <- length(GWAS$SNP)
  log_$CommonSNPs.QTL.GWAS[1] <- length(QTLs$SNP[QTLs$SNP %in% GWAS$SNP])
  
  write.csv(log_,file=paste0(out.pref,".csv"),row.names=F)
  if(flag){
    save(result,file = paste0(out.pref,".rdat"))
  }else{
    warning("There is no shared SNPs between QTLs and all tested regions in GWAS!",call. = F)
  }
sink()

}else{
  warning("There is no shared SNPs between QTLs and GWAS!",call. = F)
}
