#.....................................

#  causal inference test

#.....................................
args<-commandArgs(trailingOnly = TRUE)
# Input parameters

cpg_snp_pairs <- args[1]
genotype_data <- args[2]
methylation_data <- args[3]
phenotype_file <- args[4]
trait <- args[5]
is_binary_trait <- trimws(args[6])
covariates_num <- args[7]
covariates_fact <- args[8]
n_permutation <- as.numeric(args[9])
out_prefix <- args[10]


print("Input arguments:")
print(paste("    CpG-SNP pairs: ",cpg_snp_pairs))
print(paste("    Genotype data: ",genotype_data))
print(paste("    Methylation data: ",methylation_data))
print(paste("    Phenotype data: ",phenotype_file))
print(paste("    Trait variable: ",trait))
print(paste("    Is trait binary: ",is_binary_trait))
print(paste("    Numeric covariates: ",covariates_num))
print(paste("    Factor covariates: ",covariates_fact))
print(paste("    Number of permutation: ",n_permutation))
print(paste("    Output prefix: ",out_prefix))
##############################
library(cit)
library(stringr)
is_binary_trait <- tolower(is_binary_trait)
print("Reading and preparing inputs...")

betas <- readRDS(methylation_data)
pheno <- read.csv(phenotype_file , stringsAsFactors = F , row.names = 1)
genotype <- read.table(genotype_data , stringsAsFactors = F,row.names = 2)
colnames(genotype) <- genotype[1,]
genotype <- genotype[-1,]
cpg_snp_pairs <- read.csv(file = cpg_snp_pairs , stringsAsFactors = F )

if(!all(identical(colnames(betas),rownames(pheno)),identical(rownames(genotype),rownames(pheno)))){
  print("Rownames in phenotype and genotype data and colnames in methylation data are not matched. Making them matched...")
  index <- intersect(colnames(betas) , intersect(rownames(pheno) , rownames(genotype)))
  if(length(index)==0){
    stop("There is no common sample name in Rownames in phenotype and genotype data and colnames in methylation data")
  }else{
    betas <- betas[,index]
    genotype <- genotype[index,]
    pheno <- pheno[index,]
  }
}
colnames(genotype) <- str_remove(colnames(genotype) , pattern = "_.*")
covariates <- matrix(data = NA , nrow = 1 , ncol = 1)
if((covariates_num!="NA")&(covariates_fact!="NA")){
  
  covariates_num <- str_split_1(covariates_num,pattern = ",")
  covariates_num <- pheno[,covariates_num]
  covariates_fact <- str_split_1(covariates_fact,pattern = ",")
  covariates_fact <- pheno[,covariates_fact]
  covariates_fact <- apply(covariates_fact, 2, function(x){return(as.numeric(as.factor(x)))})
  covariates <- cbind(covariates_num, covariates_fact)
}else{
  if((covariates_num!="NA")){
    
    covariates_num <- str_split_1(covariates_num,pattern = ",")
    covariates <- pheno[,covariates_num]
  }else{
  if((covariates_fact!="NA")){
    
    covariates_fact <- str_split_1(covariates_fact,pattern = ",")
    covariates_fact <- pheno[,covariates_fact]
    covariates <- apply(covariates_fact, 2, function(x){return(as.numeric(as.factor(x)))})
  }
  }
}

if(is_binary_trait=="yes"){
  trait <- as.numeric(as.factor(pheno[,trait]))
}else{
  trait <- as.numeric(pheno[,trait])
}

print("Running CIT test...")

cpg_snp_pairs$p_cit <- NA
cpg_snp_pairs$p_TassocL <- NA
cpg_snp_pairs$p_TassocGgvnL <- NA
cpg_snp_pairs$p_GassocLgvnT <- NA
cpg_snp_pairs$p_LindTgvnG <- NA
results <- vector(mode = "list",length = nrow(cpg_snp_pairs))
names(results) <- paste0(cpg_snp_pairs[,1],"_",cpg_snp_pairs[,2])
for(i in 1:nrow(cpg_snp_pairs)){
  
  cpg <- cpg_snp_pairs[i,1]
  snp <- cpg_snp_pairs[i,2]
  
  print(paste("CpG:",cpg,",","SNP:",snp))
  
  if(is.na(covariates[1,1])){
    if(is_binary_trait)
      result = cit.bp(L = as.numeric(genotype[,snp]), G = as.numeric(betas[cpg,]), T = trait, n.perm = n_permutation)
    else
      result = cit.cp(L = as.numeric(genotype[,snp]), G = as.numeric(betas[cpg,]), T = trait, n.perm = n_permutation)
  }else{
    if(is_binary_trait)
      result = cit.bp(L = as.numeric(genotype[,snp]), G = as.numeric(betas[cpg,]), T = trait, n.perm = n_permutation)
    else
      result = cit.cp(L = as.numeric(genotype[,snp]), G = as.numeric(betas[cpg,]), T = trait, C = covariates, n.perm = n_permutation)
  }
  results[[paste0(cpg,"_",snp)]] <- result
}

print("Writing results...")

save(results, file = paste0(out_prefix,".cit.rdat"))
