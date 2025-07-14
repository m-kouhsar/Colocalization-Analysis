
Coloc.FormatQTLs <- function(QTL.file , Allel.freq.file , N.sample){
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  
  print("Reading QTL file...")
  QTLs <- fread(file=QTL.file, stringsAsFactors=F,data.table=F,header=T,nThread=getDTthreads())
  
  ## Preparing dataset for coloc analysis
  ## Ref: https://stats.stackexchange.com/questions/337070/compute-standard-error-from-beta-p-value-sample-size-and-the-number-of-regres
  
  QTLs$SE <- QTLs$beta/abs(QTLs$statistic)
  QTLs$SNPID <- paste(str_split(QTLs$snps,pattern=":",simplify=T)[,1],str_split(QTLs$snps,pattern=":",simplify=T)[,2],sep = ":")
  QTLs$SNP.Pos <- str_split(QTLs$SNPID,pattern=":",simplify=T)[,2]
  QTLs$CHR <- str_split(QTLs$SNPID,pattern=":",simplify=T)[,1]
  
  print("Reading allele frequency file...")
  maf <- fread(file=Allel.freq.file, stringsAsFactors=F,data.table=F,header=T,nThread=getDTthreads())
  
  print("Adding MAF to QTLs...")
  snps <- str_split(QTLs$snps, pattern="_" , simplify=T)[,1]
  index <- match(snps , maf$SNP)
  QTLs$MAF <- maf$MAF[index]
  
  print("Adding number of samples column...")
  QTLs$N <- N.sample
  head(QTLs)
  
  print("Saving results...")
  QTL.file <- str_remove(QTL.file,pattern=".csv")
  write.csv(QTLs,file=paste0(QTL.file,".formatted.csv"),row.names=F)
}
##GWAS$beta <- GWAS$Effect <- log(OR) #beta = ln(odd ratio)
##GWAS$maf <- ((GWAS$fcas*GWAS$ncas)+(GWAS$fcon*GWAS$ncon))/(GWAS$ncas+GWAS$ncon)
