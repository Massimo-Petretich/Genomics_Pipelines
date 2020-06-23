#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/featureCounts_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/featureCounts_err.txt



bamFiles=''

outputDirectory=''
tagOutput='counts'

cores=4

annotation='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Annotations/hg19_Ensembl_GRCh37.gtf'
strandSpecific=0	# 0 (unstranded), 1 (stranded) for example Lexogen forward, 2 (reversely stranded) for example for the NEX next Ultra II stranded kit; for the featureCounts function
GTFattrType='gene_id'	# default for transcriptome assembly; for the featureCounts function

lib='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343'



if [ ! -d $outputDirectory ];                   then    mkdir $outputDirectory;                 fi

source activate r-343

R     --slave     --args     $annotation     $bamFiles       $outputDirectory/$tagOutput.RData	$cores	$strandSpecific	$GTFattrType	$lib  <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/featureCounts.R


exit
