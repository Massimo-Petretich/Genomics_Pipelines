#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/Table_Count_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/Table_Count_err.txt


files=''		# .txt file with 1 file path per line (Granges RData files)
GRangesObject=''	# .RData file containing the query regions
outputFile=''		# .RData object (GRanges)
cores=4

source activate r-343

R     --slave     --args     $GRangesObject $files	$outputFile	$cores  <	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/table_count.R

