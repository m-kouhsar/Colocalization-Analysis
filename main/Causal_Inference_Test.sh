#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

### print start date and time
echo Job started on:
date -u

set -e
####### 

### NOTE: Do not store confidenial information in this file use the config file

######
SCRIPTDIR=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/scripts
cpg_snp_pairs=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/inputs/cit/SNPs.CpGs.csv
genotype_data=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/inputs/cit/EMIF_AD_imputed.ClumpedSNPs.raw
methylation_data=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/inputs/cit/EMIF_AD_Beta.SelectedCpGs.rds
phenotype_file=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/inputs/cit/EMIF_AD_Pheno.csv
trait=Central_CSF_YKL40
is_binary_trait=0            #1 means True
covariates_num=NA
covariates_fact=NA
n_permutation=500
out_prefix=/lustre/projects/Research_Project-191391/Morteza/coloc.Ehsan/results/EMIF_AD


module load R

Rscript ${SCRIPTDIR}/Causal_Inference_Test.R $cpg_snp_pairs $genotype_data $methylation_data $phenotype_file $trait $is_binary_trait $covariates_num $covariates_fact $n_permutation $out_prefix


echo Job finished:
date -u
