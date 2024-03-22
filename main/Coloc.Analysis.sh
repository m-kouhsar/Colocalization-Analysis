#!/bin/sh
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-11

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Scripts/Coloc

gwas_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/summary_stat
loci_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/summary_stat/regions
qtl_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc

out_prefix=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc/darkgreen/darkgreen.cis

ref_genome_prefix=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Ref/g1000_eur/g1000_eur_rsid

# Ref genome will use to calculate LD for the lead SNPs in loci file. So, the SNP IDs in Ref genome must be the same as lead SNPs

input_list=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/March2024/Results/Coloc/coloc.input.txt

distance_=1000 #kb
use_ld=yes
ld_threshold=0.6
type_="cc" # quant for quantitative trait and cc for binary

###################################################################################
mkdir -p "$(dirname "${out_prefix}")"

declare -A dupp_loci

while read -r line
do
  IFS=',' read -r -a array <<< "$line"
  gwas+=(${array[0]})
  loci+=(${array[1]})
  qtl+=(${array[2]})
  ((dupp_loci[$loci]++))
done < $input_list

if [ $use_ld = "yes" ]
then
  out_prefix="${out_prefix}.ColocResult.${type_}.dist.${distance_}.ld.${ld_threshold}.Output$(( SLURM_ARRAY_TASK_ID + 1 ))"
else
  out_prefix="${out_prefix}.ColocResult.${type_}.dist.${distance_}.Output$(( SLURM_ARRAY_TASK_ID + 1 ))"
fi

loci_file=${loci_dir}/${loci[$SLURM_ARRAY_TASK_ID]}
gwas_file=${gwas_dir}/${gwas[$SLURM_ARRAY_TASK_ID]}
qtl_file=${qtl_dir}/${qtl[$SLURM_ARRAY_TASK_ID]}

#if [  ${dupp_loci[${loci[$SLURM_ARRAY_TASK_ID]}]} -gt 1 ]
#then
#	cp ${loci_dir}/${loci[$SLURM_ARRAY_TASK_ID]} ${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
#  loci_file=${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
#fi

Rscript ${ScriptDir}/Coloc.Analysis.V1.R   $qtl_file  $gwas_file $loci_file $type_ $distance_ $use_ld $ld_threshold $ref_genome_prefix $out_prefix  

#if [  ${dupp_loci[${loci[$SLURM_ARRAY_TASK_ID]}]} -gt 1 ]
#then
#	rm ${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
#fi

echo "Done!"
