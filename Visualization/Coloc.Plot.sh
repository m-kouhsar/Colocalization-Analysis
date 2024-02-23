#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=4:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


qtl_file=./QTL.ColocPlot.csv
cpg_file=./CpGs.ColocPlot.txt
gwas_file=./AD.kunkle.ColocPlot.txt
gwas_trait=Alzheimer’s disease
out_prefix=./Coloc.Plot
cpg_interest=cg07050504,cg24774200,cg15831875
pvalue_GWAS=1e-5
pvalue_QTL=1e-5
range=1e+6
gbuild=hg19
congruence=F

ScriptDir=./Visualization

Rscript ${ScriptDir}/Plot.Coloc.R $qtl_file $cpg_file $gwas_file "$gwas_trait" $out_prefix $cpg_interest $pvalue_GWAS $pvalue_QTL $range $gbuild $congruence
