#!/bin/sh
#SBATCH -A Research_Project-MRC164847 
#SBATCH --export=ALL 
#SBATCH -D . 
#SBATCH -p mrcq
#SBATCH --time=48:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mail-type=END 
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk 

########################################################################################################################################
######################################### Arguments description ########################################################################

#  input_list: A csv file with 4 columns (1st: GWAS Summary statistics file name, 2nd: loci file name, 3rd: QTL file name, 4th: Output file prefix). 
#              A coloc analysis will be performed for each row in this file. 
#  gwas_dir: A directory contains the GWAS summary statistic files (first column in the input_list)
#            The GWAS summary statistic file must be a matrix like file contains the following columns (column names are not case sensitive):
#            SNP: SNP IDs. Must be matched with the SNP IDs in QTLs
#            CHR: Chromosome
#            POS: Position
#            MAF: Minor allele frequency
#            SE: Standard error
#            BETA: Effect size
#            N: Number of samples
#  loci_dir: A directory contains the GWAS locus/regions informations (second column in the input_list).
#            locus information file is a matrix like file contains the following columns (column names are not case sensitive):
#            ID: An ID for the region of interest (this is an optional column. An ID will be generated if it doesn't specified)
#            CHR: Region chromosome
#            START: Region start position
#            END: Region end position
#  qtl_dir: A directory contains the QTLs informations (3rd column in the input_list).
#           The QTL file is a matrix like file contains the following columns (column names are not case sensitive):
#           SNP: SNP IDs. Must be matched with the SNP IDs in the GWAS file
#           Gene: The gene/CpG IDs
#           MAF: SNPs minor allele frequency
#           SE: Standard error
#           BETA: Effect size
#           N: Number of samples
#  use_ld: In the locus/region information file, you can specify the lead SNPs information instead of regions (ID, CHR, START, END). 
#          So, if you set this parameter to "yes", for all regions with length < 2 (including lead SNPs), the linkage disequilibrium will be calculated 
#          based on a reference genome genotype data and then a region will be defined based on a LD threshold and a distance (window size). 
#  ref_genome_prefix: plink binary prefix for the reference genome. A refrence genome genotype data in binary plink format is requiered if you set use_ld to "yes"
#  ld_threshold: LD threshold to define the regions based on LD
#  distance_: A distance (in KB) to extend the regions with length < 2
#             If you set use_ld to "yes", this distance will be used as "--ld-window-kb" option in "plink --r2" command
#             If you set use_ld to "no", this distance will be simply use to define a region aroud the lead SNPs [START - distance_ , END + distance_]
#  type_QTL: Type of trait in QTL data ('quant' for quantitative trait and 'cc' for binary). For more infromation see the coloc package manual. 
#  type_GWAS: Type of trait in GWAS data ('quant' for quantitative trait and 'cc' for binary). For more infromation see the coloc package manual. 
#  use_permut: If you set this option to "yes", A permutation test will be perform on each Significant Coloc Results (SCR). 
#              We have an analysis for each QTL-Region pair. So, for each SCR QTL-Region pair, 1000 random region with the same length (or the same way for regions with length < 2) 
#              will be generated. After getting the results on all of these random regions, a P-value will be calculated using the following formula:
#              Pvalue = (number of random SCR)/1000. 
#              A Significant Coloc Result (SCR) is a results contains at least one probabilty (H0-H4 or H3+H4) larger than coloc_P_threshold. 
#  coloc_P_threshold=0.9 A probabilty threshold to define Significant Coloc Results
#  n_permut=1000
#  out_prefix: Output files prefix
#  ScriptDir: A directory contains all dependant scripts

###########################################################################################################################################

input_list=./coloc.input.csv 
gwas_dir=./summary_stat
loci_dir=./regions
qtl_dir=./QTLs
use_ld=yes
ref_genome_prefix=./g1000_eur_ChrPos
ld_threshold=0.6
distance_=250
type_QTL="cc"
type_GWAS="cc"
use_permut=yes 
coloc_P_threshold=0.9
n_permut=1000
out_prefix=./Results/Sample_coloc
ScriptDir=./Main/Dependant_Scripts

###########################################################################################################################################

dos2unix $input_list
line_count=$(wc -l < $input_list)

sbatch --array=0-$(($line_count - 1)) ${ScriptDir}/Coloc.Permutation.sh $input_list $gwas_dir $loci_dir $qtl_dir $use_ld $ref_genome_prefix $ld_threshold $distance_ $type_QTL $type_GWAS $use_permut $coloc_P_threshold $n_permut $out_prefix $ScriptDir
