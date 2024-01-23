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
#SBATCH --array=0-2  ## Number of lines in coloc.inputs.txt file

module load R

ScriptDir=./main

SumStat_dir=./SamStats
loci_dir=./Regions
qtl_dir=./QTLs

out_dir=./Outputs

ref_genome_prefix=./Ref/g1000_eur/g1000_eur_rsid

input_list=./coloc.inputs.txt

ld_threshold=0.6
distance_=250000 #bp
type_="quant"

###################################################################################
mkdir -p $out_dir

declare -A dupp 

while read -r line
do
	IFS=',' read -r -a array <<< "$line"
	gwas+=(${array[0]})
	qtl+=(${array[1]})
	loci+=(${array[2]})
	((dupp[$loci]++))
done < $input_list

out_pref1=${gwas[$SLURM_ARRAY_TASK_ID]%".txt"}
out_pref2=${qtl[$SLURM_ARRAY_TASK_ID]%".csv"}
out_pref=${out_pref1}_${out_pref2}

loci_file=${loci_dir}/${loci[$SLURM_ARRAY_TASK_ID]}

if [  ${dupp[${loci[$SLURM_ARRAY_TASK_ID]}]} -gt 1 ]
then
	cp ${loci_dir}/${loci[$SLURM_ARRAY_TASK_ID]} ${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
  loci_file=${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
fi

Rscript ${ScriptDir}/Coloc.Analysis.R   ${qtl_dir}/${qtl[$SLURM_ARRAY_TASK_ID]}  ${SumStat_dir}/${gwas[$SLURM_ARRAY_TASK_ID]} $loci_file $ref_genome_prefix $distance_ $type_ $ld_threshold ${out_dir}/$out_pref

if [  ${dupp[${loci[$SLURM_ARRAY_TASK_ID]}]} -gt 1 ]
then
	rm ${out_dir}/${loci[$SLURM_ARRAY_TASK_ID]%".csv"}.${SLURM_ARRAY_TASK_ID}".csv"
fi

echo "Done!"
