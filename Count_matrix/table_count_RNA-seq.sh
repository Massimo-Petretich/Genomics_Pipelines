#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/Table_Count_RNA-seq_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/Table_Count_RNA-seq_err.txt


files=''		# .txt file cintaining per line 1 path to a .bam files
outputFile=''		# .RData output file (table of counts)

cores=4

ignoreStrand='TRUE'

singleEnd='TRUE'

mode="Union"


annotation='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Annotations/hg19_UCSC.RData'


source activate r-343
R     --slave     --args     $annotation $files	$outputFile	$cores  $ignoreStrand	$singleEnd	$mode	<	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/table_count_RNA-seq.R
