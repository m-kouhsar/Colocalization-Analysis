#!/bin/sh
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-11  # Numer of array should be equal to number of lines in coloc.input.csv file (each line will be run by an array)

ScriptDir=./Scripts/Coloc

gwas_dir=./summary_stat
loci_dir=./regions
qtl_dir=./QTLs

out_prefix=./Results.Coloc # csv file with 4 columns (1st: GWAS Summary statistics file name, 2nd: loci file name, 3rd: QTL file name, 4th: Output file ID)

ref_genome_prefix=./g1000_eur_rsid

# Ref genome will use to calculate LD for the lead SNPs in loci file. So, the SNP IDs in Ref genome must be the same as lead SNPs

input_list=./coloc.input.csv

distance_=1000 #kb
use_ld=yes
ld_threshold=0.6
type_="cc" # quant for quantitative trait and cc for binary

###################################################################################
# conda install conda-forge::dos2unix
dos2unix $input_list
mkdir -p "$(dirname "${out_prefix}")"

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

Rscript ${ScriptDir}/Coloc.Analysis.V1.R   $qtl_file  $gwas_file $loci_file $type_ $distance_ $use_ld $ld_threshold $ref_genome_prefix $out_prefix  

echo "Done!"
