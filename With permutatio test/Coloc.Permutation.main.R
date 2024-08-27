
args<-commandArgs(TRUE)

QTL.file <- args[1]
GWAS.file <- args[2]
loci.file <- args[3]
type_QTL <- args[4]                         ## "quant" for quantitative traits and "cc" for binary traits
type_GWAS <- args[4]                         ## "quant" for quantitative traits and "cc" for binary traits
distance <- as.numeric(args[5])          ## distance in bp to define regions based on LD (it will be ignored if length locus > 2)
LD <- ifelse(tolower(args[6])=="yes",T,F)
LD.threshold <- as.numeric(args[7])      ## LD threshold for defining regions (it should be between 0 and 1)
ref.genome.prefix <- args[8]
use.permut <- ifelse(tolower(args[9])=="yes",T,F)
n.permut <- as.numeric(args[10])
p.threshold <- as.numeric(args[11])
out.pref <- args[12]
log <- file(paste0(out.pref , ".log.txt") , open = "wt")
sink(file = log , type = "output")
sink(file = log , type = "message")
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
cat("     QTL trait type= ",type_QTL,"\n")
cat("     GWAS trait type= ",type_GWAS,"\n")
cat("     Distance for selecting region based on LD (if the region length < 2)= ",distance, " KB\n")
cat("     Using LD to define regions for lead SNPs? ",ifelse(LD,"Yes","No"),"\n")
cat("     LD threshold= ",LD.threshold,"\n")
cat('\n')
cat("     Reference genome binary files name (--bfile option in plink)= ",ref.genome.prefix,"\n")
cat("     Output prefix= ",out.pref,"\n")
cat('\n')
cat("     Using random regions for calculating result P-value? ",ifelse(use.permut,"Yes","No"),"\n")
cat("    Coloc probability threshold= ",p.threshold,"\n")
cat("     Number of random regions= ",n.permut,"\n")
cat("#############################################################################################\n")
cat("\n")
cat("Loading libraries...\n")
suppressMessages(library(coloc))
suppressMessages(library(stringr))
suppressMessages(library(vroom))
suppressMessages(library(data.table))
options(datatable.fread.datatable=FALSE)
suppressMessages(library("funr"))

source(paste0(dirname(sys.script()),"/Coloc.Permutation.Functions.R"))

cat("Reading QTL file...\n")
QTLs <- vroom(QTL.file, col_types = c(.default = "c"),num_threads=16)
QTLs <- as.data.frame(QTLs)
names(QTLs) <- toupper(names(QTLs))

QTLs.col <- c("SNP","MAF","SE","BETA","N","GENE")
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

GWAS.col <- c("SNP","MAF","SE","BETA","N","CHR","POS")
if(!all(GWAS.col %in% colnames(GWAS))){
  missed.col <- GWAS.col[which(!(GWAS.col %in% colnames(GWAS)))]
  stop("The following columns are missed in GWAS data: ",paste(missed.col , collapse = ","))
}

cat("Removing dupplicated SNPs from GWAS summary statistics...\n")
GWAS <- GWAS[!duplicated(GWAS$SNP),]
GWAS$MAF <- as.numeric(GWAS$MAF)
GWAS$SE <- as.numeric(GWAS$SE)
GWAS$BETA <- as.numeric(GWAS$BETA)
GWAS$N <- as.numeric(GWAS$N)
GWAS$POS <- as.numeric(GWAS$POS)
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
  colnames(loci) <- toupper(colnames(loci))
  loci.col <- c("CHR","START","END")
  if(!all(loci.col %in% colnames(loci))){
    missed.col <- loci.col[which(!(loci.col %in% colnames(loci)))]
    stop("The following columns are missed in GWAS data: ",paste(missed.col , collapse = ","))
  }
  loci$START <- as.numeric(loci$START)
  loci$END <- as.numeric(loci$END)
  if(!("ID" %in% names(loci))){
    if("RSID" %in% names(loci)){
      loci$ID <- loci$RSID
    }else{
      if("SNP" %in% names(loci)){
        loci$ID <- loci$SNP
      }else{
        loci$ID <- paste0("loci.",c(1:nrow(loci)))
      }
    }
  }
  
  need.definition <- (loci$END - loci$START) < 2
  if(any(need.definition)){
    cat(sum(need.definition)," lead SNPS (loci with lenghth < 2) were detected.\n")
    cat("Generating regions for lead SNPs...\n")
    if(!("SNP" %in% names(loci))){
       if(!("RSID" %in% names(loci))){
         stop("SNP or RSID not found in loci file columns!")
       }else{
         loci$SNP <- loci$RSID
       }
    }
    if(LD){
      cat("Calculatink LD using plink...\n")
      ld.data <- calculate.LD(ld.window.r2 = LD.threshold , ld.window.kb = distance,
                              snps.list = loci$SNP[need.definition],ref.genome.prefix = ref.genome.prefix,out.prefix = out.pref)
      for (i in 1:sum(need.definition)) {
        if(loci$SNP[i] %in% ld.data$SNP_A){
          ld1 <- ld.data[ld.data$SNP_A==loci$SNP[i],]
          loci$START[i] <- min(ld1$BP_B)
          loci$END[i] <- max(ld1$BP_B)
        }
      }
    }else{
      loci$START[need.definition] <- loci$START[need.definition] - (distance*1000)
      loci$END[need.definition] <- loci$END[need.definition] + (distance*1000)
    }
  }
  
  log_ <- as.data.frame(matrix(data = NA,nrow = nrow(loci),ncol = 30))
  names(log_) <- c("Coloc.Type","LD.threshold","Distance","UniqGene.QTL","UniqSNPs.QTL","UniqSNPs.GWAS","CommonSNPs.QTL.GWAS",
                   "locus.Id","locus.Chr","locus.Start","locus.End","locus.length","UniqSNPs.locus","CommonSNPs.QTL.locus","PP.H0.abf","PP.H1.abf",
                   "PP.H2.abf","PP.H3.abf","PP.H4.abf","sum.H3.H4","SNPs","nSNPs","Genes","nGenes","H0.Pvalue","H1.Pvalue","H2.Pvalue","H3.Pvalue","H4.Pvalue","Sum.H3.H4.Pvalue")
  
  result <- vector(mode = "list",length = nrow(loci))
  names(result) <- paste0("locus.",loci$ID)
  permut.results <- list()
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
    log_$locus.Chr[i] = loci$CHR[i]
    log_$locus.Start[i] <- loci$START[i]
    log_$locus.End[i] <- loci$END[i]
    log_$locus.length[i] <- loci$END[i] - loci$START[i]
    index <- (GWAS$CHR == loci$CHR[i]) & (GWAS$POS > loci$START[i]) & (GWAS$POS < loci$END[i])
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
      
      if(use.permut & any(log_$PP.H0.abf[i] > p.threshold, log_$PP.H1.abf[i] > p.threshold, log_$PP.H2.abf[i] > p.threshold, 
                          log_$PP.H3.abf[i] > p.threshold, log_$PP.H4.abf[i] > p.threshold, log_$sum.H3.H4[i] > p.threshold)){
        cat("Generating random regions...\n")
        loci.permut <- permutation(locus.snp = loci$SNP[i],locus.start = loci$START[i],locus.end = loci$END[i], LD = LD,
                                   LD.threshold = LD.threshold,
                                   distance = distance,ref.genome.prefix = ref.genome.prefix, out.prefix = out.pref,
                                   snps.all = subset(GWAS,select = c("SNP","CHR","POS")),n = n.permut)
        loci.permut$PP.H0.abf <- 0
        loci.permut$PP.H1.abf <- 0
        loci.permut$PP.H2.abf <- 0
        loci.permut$PP.H3.abf <- 0
        loci.permut$PP.H4.abf <- 0
        loci.permut$sum.H3.H4 <- 0
        cat("Running coloc on random regions...\n")
        for (j in 1:nrow(loci.permut)) {
          cat("Working on random regions ",j,"...\n")
          index <- (GWAS$CHR == loci.permut$CHR[j]) & (GWAS$POS > loci.permut$START[j]) & (GWAS$POS < loci.permut$END[i])
          GWAS.loci <- GWAS[index , ]
          c <- length(intersect(QTLs$SNP , GWAS.loci$SNP))
          if(c > 0){
            GWAS_coloc = list(beta = GWAS.loci$BETA,
                              varbeta = GWAS.loci$SE^2,
                              type = type_,
                              snp = GWAS.loci$SNP, 
                              MAF = GWAS.loci$MAF,
                              N = GWAS.loci$N)
            result.coloc = coloc.abf(GWAS_coloc, QTLs_coloc)
            temp = as.data.frame(result.coloc$summary)
            loci.permut$PP.H0.abf[j] = temp["PP.H0.abf",1]
            loci.permut$PP.H1.abf[j] = temp["PP.H1.abf",1]
            loci.permut$PP.H2.abf[j] = temp["PP.H2.abf",1]
            loci.permut$PP.H3.abf[j] = temp["PP.H3.abf",1]
            loci.permut$PP.H4.abf[j] = temp["PP.H4.abf",1]
            loci.permut$sum.H3.H4[j] = temp["PP.H3.abf",1] + temp["PP.H4.abf",1]
          }
        }
        permut.results <- append(permut.results,list(loci.permut))
        names(permut.results)[length(permut.results)] <- loci$ID[i]
        
        log_$H0.Pvalue[i] <- 1-sum(p.threshold > loci.permut$PP.H0.abf )/nrow(loci.permut)
        log_$H1.Pvalue[i] <- 1-sum(p.threshold > loci.permut$PP.H1.abf )/nrow(loci.permut)
        log_$H2.Pvalue[i] <- 1-sum(p.threshold > loci.permut$PP.H2.abf )/nrow(loci.permut)
        log_$H3.Pvalue[i] <- 1-sum(p.threshold > loci.permut$PP.H3.abf )/nrow(loci.permut)
        log_$H4.Pvalue[i] <- 1-sum(p.threshold > loci.permut$PP.H4.abf )/nrow(loci.permut)
        log_$Sum.H3.H4.Pvalue[i] <- 1-sum(p.threshold > loci.permut$sum.H3.H4 )/nrow(loci.permut)
        cat("H0 P-value: ",log_$H0.Pvalue[i],"\n")
        cat("H1 P-value: ",log_$H1.Pvalue[i],"\n")
        cat("H2 P-value: ",log_$H2.Pvalue[i],"\n")
        cat("H3 P-value: ",log_$H3.Pvalue[i],"\n")
        cat("H4 P-value: ",log_$H4.Pvalue[i],"\n")
        cat("(H3 + H4) P-value: ",log_$Sum.H3.H4.Pvalue[i],"\n")
        cat("\n")
      }
    }
    cat('\n')
  }
  
  log_ <- log_[order(log_$sum.H3.H4 , decreasing = T),]
  log_$Coloc.Type[1] <- type_
  log_$LD.threshold[1] <- LD.threshold
  log_$Distance[1] <- distance
  log_$UniqGene.QTL[1] <- length(unique(QTLs$GENE))
  log_$UniqSNPs.QTL[1] <- length(QTLs$SNP)
  log_$UniqSNPs.GWAS[1] <- length(GWAS$SNP)
  log_$CommonSNPs.QTL.GWAS[1] <- length(QTLs$SNP[QTLs$SNP %in% GWAS$SNP])
  
  write.csv(log_,file=paste0(out.pref,".csv"),row.names=F)
  if(flag){
    save(result,file = paste0(out.pref,".rdat"))
    save(permut.results,file = paste0(out.pref,".permutation.rdat"))
  }else{
    warning("There is no shared SNPs between QTLs and all tested regions in GWAS!",call. = F)
  }
  
}else{
  warning("There is no shared SNPs between QTLs and GWAS!",call. = F)
}
sink()
#sink()
