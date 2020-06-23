#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/fragmentsFinder_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/fragmentsFinder_err.txt


region='' # UCSC browser formatted genomic coordinate
pattern_1='' # restriction enzyme palindromic sequence (1st cut like 'GATC')
pattern_2='' # restriction enzyme palindromic sequence (1st cut like 'GTAC')
output_directory=''
tag_output=''

cores=8

BSgenome_package='BSgenome.Hsapiens.UCSC.hg19'  

lib='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343'



if [ ! -d $output_directory ];                   then    mkdir $output_directory;                 fi

if [ ! -d $output_directory/$tag_output ];       then    mkdir $output_directory/$tag_output;     fi

output_directory=$output_directory/$tag_output
cd $output_directory


source activate r-343


R  --slave --args $region $pattern_1 $pattern_2 $BSgenome_package $output_directory $cores $lib < /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/fragmentsFinder.R

exit
