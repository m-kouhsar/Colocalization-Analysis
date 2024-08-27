#!/bin/sh
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=48:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-1 ## Numer of array should be equal to number of lines in coloc.input.csv file (each line will be run by an array)


ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Scripts/Coloc

gwas_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/summary_stat
loci_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/summary_stat/regions
qtl_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc

out_prefix=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc/darkgreen/darkgreen.Cis

ref_genome_prefix=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Ref/g1000_eur/g1000_eur_ChrPos

# Ref genome will use to calculate LD for the lead SNPs in loci file. So, the SNP IDs in Ref genome must be the same as lead SNPs

input_list=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc/coloc.input.csv # csv file with 4 columns (1st: GWAS Summary statistics file name, 2nd: loci file name, 3rd: QTL file name, 4th: Output file ID)


distance_=250 #kb
use_ld=yes
ld_threshold=0.6
coloc_P_threshold=0.9 #To Save time, for any results (H0, H1, H2,H3,H4, and H3+H4) larger than this threshold the permutation test will run
type_QTL="cc" # 'quant' for quantitative trait and 'cc' for binary
type_GWAS="cc" # 'quant' for quantitative trait and 'cc' for binary
use_permut=yes 
n_permut=1000

###################################################################################
dos2unix $input_list
mkdir -p "$(dirname "${out_prefix}")"


while read -r line
do
	IFS=',' read -r -a array <<< "$line"
	gwas+=(${array[0]})
	loci+=(${array[1]})
	qtl+=(${array[2]})
   OutID+=(${array[3]})

done < $input_list

if [ $use_ld = "yes" ]
then
  out_prefix="${out_prefix}.${OutID[$SLURM_ARRAY_TASK_ID]}.QTL.${type_QTL}.GWAS.${type_GWAS}.dist.${distance_}.ld.${ld_threshold}.Coloc"
else
  out_prefix="${out_prefix}.${OutID[$SLURM_ARRAY_TASK_ID]}.QTL.${type_QTL}.GWAS.${type_GWAS}.dist.${distance_}.Coloc"
fi

if [ $use_permut = "yes" ]
then
  out_prefix="${out_prefix}.Permutation.${n_permut}.Output"
else
  out_prefix="${out_prefix}.Output"
fi

loci_file=${loci_dir}/${loci[$SLURM_ARRAY_TASK_ID]}
gwas_file=${gwas_dir}/${gwas[$SLURM_ARRAY_TASK_ID]}
qtl_file=${qtl_dir}/${qtl[$SLURM_ARRAY_TASK_ID]}



Rscript ${ScriptDir}/Coloc.Permutation.main.R   $qtl_file  $gwas_file $loci_file $type_QTL $type_GWAS $distance_ $use_ld $ld_threshold $ref_genome_prefix $use_permut $n_permut $coloc_P_threshold $out_prefix  


echo "Done!"
