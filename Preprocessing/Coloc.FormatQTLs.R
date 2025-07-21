
Coloc.FormatQTLs <- function(QTLs , Allel.freq , N.sample){
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  
  if(!(is.data.frame(QTLs) | is.matrix(QTLs))){
    print("Reading QTL file...")
    QTLs <- fread(file=QTLs, stringsAsFactors=F,data.table=F,header=T,nThread=getDTthreads())
  }
  
  if(!(is.data.frame(Allel.freq) | is.matrix(Allel.freq))){
    print("Reading allele frequency file...")
    Allel.freq <- fread(file=Allel.freq, stringsAsFactors=F,data.table=F,header=T,nThread=getDTthreads())
  }
  
  ## Preparing dataset for coloc analysis
  ## Ref: https://stats.stackexchange.com/questions/337070/compute-standard-error-from-beta-p-value-sample-size-and-the-number-of-regres
  message("Converst SNP IDs...")
  QTLs$SE <- QTLs$beta/abs(QTLs$statistic)
  QTLs$SNP <- paste(str_split(QTLs$snps,pattern=":",simplify=T)[,1],str_split(QTLs$snps,pattern=":",simplify=T)[,2],sep = ":")
  QTLs$SNP.Pos <- str_split(QTLs$SNPID,pattern=":",simplify=T)[,2]
  QTLs$CHR <- str_split(QTLs$SNPID,pattern=":",simplify=T)[,1]
  
  message("Adding MAF to QTLs...")
  snps <- str_split(QTLs$snps, pattern="_" , simplify=T)[,1]
  index <- match(snps , Allel.freq$SNP)
  QTLs$MAF <- Allel.freq$MAF[index]
  
  message("Adding number of samples column...")
  QTLs$N <- N.sample
  
  return(QTLs)
}
##GWAS$beta <- GWAS$Effect <- log(OR) #beta = ln(odd ratio)
##GWAS$maf <- ((GWAS$fcas*GWAS$ncas)+(GWAS$fcon*GWAS$ncon))/(GWAS$ncas+GWAS$ncon)
