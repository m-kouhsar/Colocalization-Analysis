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
######
SCRIPTDIR=./main
cpg_snp_pairs=./SNPs.CpGs.csv
genotype_data=./Genotype.raw
methylation_data=./Methylation.Beta.Values.rds
phenotype_file=./Pheno.csv
trait=AD
is_binary_trait=1            #1 means True
covariates_num=Age,BMI
covariates_fact=Sex,BraakStage
n_permutation=500
out_prefix=./Results.CIT


module load R

Rscript ${SCRIPTDIR}/Causal_Inference_Test.R $cpg_snp_pairs $genotype_data $methylation_data $phenotype_file $trait $is_binary_trait $covariates_num $covariates_fact $n_permutation $out_prefix


echo Job finished:
date -u
