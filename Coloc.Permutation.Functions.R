
###############################################################################################
calculate.LD <- function(ld.window.r2 = 0.6 , ld.window.kb , snps.list , ref.genome.prefix,out.prefix){
  write.table(snps.list,file =paste0(out.prefix,".LociLeadSNPs.txt"),col.names = F , row.names = F,quote = F )
  plink_cmd <- paste("plink --r2 --bfile",ref.genome.prefix , "--ld-snp-list", paste0(out.prefix,".LociLeadSNPs.txt") ,"--ld-window-r2",ld.window.r2,"--ld-window-kb",ld.window.kb,"--out",out.prefix)
  exit_code <- system(plink_cmd, ignore.stdout=T,wait = T)
  if(!(exit_code>0)){
    ld <- fread(paste0(out.prefix,".ld"),stringsAsFactors = F , header = T)
    #rm <- file.remove(paste0(out.prefix,".log"))
    #rm <- file.remove(paste0(out.prefix,".ld"))
    #rm <- file.remove(paste0(out.prefix,".SNPList.txt"))
  }else{
    #rm <- file.remove(paste0(out.prefix,".SNPList.txt"))
    stop("Could not calculate LD using plin!","\n","plink command:","\n",plink_cmd,"\n","exit code:",exit_code)
  }
  return(ld)
}

####################################################################################
permutation <- function(locus.snp,locus.start,locus.end, LD = T, LD.threshold = 0.6 , distance=1000 , ref.genome.prefix, 
                        snps.all, n=1000, out.prefix){
  loci.permut <- as.data.frame(matrix(data = 0 , nrow = n , ncol = 4))
  colnames(loci.permut) <- c("SNP","CHR","START","END")
  
  snps.permut <- snps.all[sample(1:nrow(snps.all) , size = n,replace = F),]
  if((locus.end - locus.start) <2 ){
    if(LD){
      ld.data <- calculate.LD(ld.window.r2 = LD.threshold , ld.window.kb = distance,
                              snps.list = snps.permut$SNP,ref.genome.prefix = ref.genome.prefix,out.prefix = out.prefix)
    }
    for (i in 1:n) {
      loci.permut$SNP[i] <- snps.permut$SNP[i]
      loci.permut$CHR[i] <- snps.permut$CHR[i]
      loci.permut$START[i] <- snps.permut$POS[i]
      loci.permut$END[i] <- snps.permut$POS[i]
      if(LD){
        if(loci.permut$SNP[i] %in% ld.data$SNP_A){
          ld1 <- ld.data[ld.data$SNP_A==loci.permut$SNP[i],]
          loci.permut$START[i] <- min(ld1$BP_B)
          loci.permut$END[i] <- max(ld1$BP_B)
        }
      }else{
        loci.permut$START[i] <- loci$START[i] - (distance*1000)
        loci.permut$END[i] <- loci$END[i] + (distance*1000)
      }
    }
  }else{
    L <- round((locus.end - locus.start)/2,digits = 0)
    for (i in 1:n) {
      loci.permut$SNP[i] <- snps.permut$SNP[i]
      loci.permut$CHR[i] <- snps.permut$CHR[i]
      loci.permut$START[i] <- snps.permut$POS[i] - L
      loci.permut$END[i] <- snps.permut$POS[i] + L
    }
  }
  return(loci.permut)
}
