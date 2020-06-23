#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10-10:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/motif_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/motif_err.txt



peaks=''
mainDir=''		# home directory of the ChIP not the one of the peak calling
tag=''			# prefix of the file generated
percIdentity=85		# [80 | 85 | 90] for mm and  hs, standard = 85
matrices='Human_PWM_full_collection'	# filename in package data folder without .rda
cores=4
lib='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343'





source activate r-343

R     --slave     --args	$peaks	$mainDir	$tag	$percIdentity	$matrices	$cores	$lib  <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/motif2.R


exit
