#!/bin/sh
#SBATCH -A Research_Project-MRC164847 
#SBATCH --export=ALL 
#SBATCH -D . 
#SBATCH -p mrcq
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mail-type=END 
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk 

#################################################################################################################################################################################################################

echo "Job started at $(date '+%H:%M:%S') on $(date '+%d/%m/%Y')"

input_list=$1
gwas_dir=$2
loci_dir=$3
qtl_dir=$4
use_ld=$5
ref_genome_prefix=$6
ld_threshold=$7
distance_=$8
type_QTL=$9
type_GWAS=${10}
use_permut=${11}
coloc_P_threshold=${12}
n_permut=${13}
out_prefix=${14}
ScriptDir=${15}


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

echo "Find the logs in ${out_prefix}.log.txt"

Rscript ${ScriptDir}/Coloc.Permutation.R   $qtl_file  $gwas_file $loci_file $type_QTL $type_GWAS $distance_ $use_ld $ld_threshold $ref_genome_prefix $use_permut $n_permut $coloc_P_threshold $out_prefix  


echo "Job finished at $(date '+%H:%M:%S') on $(date '+%d/%m/%Y')"
