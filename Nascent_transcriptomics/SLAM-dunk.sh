#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/SLAM-dunk_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/SLAM-dunk_err.txt

cores=16

fastq_files=''
output_directory=''
output_file=''
max_ylim=5

bed_file='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Annotations/hg38_refseq_ensembl_3UTR_fillup.bed'
genome='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa'


read_length=$(gunzip -c $fastq_files | head -n 4 | awk 'FNR==2{print $1}' | wc | awk '{print $3}')


if [ ! -d $output_directory ]; then mkdir $output_directory; fi
if [ ! -d $output_directory/$output_file ];       then    mkdir $output_directory/$output_file;     fi

output_directory=$output_directory/$output_file

cd $output_directory

gunzip -c $fastq_files > $output_file


/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $output_file   --outdir=$output_directory


################ SLAM-dunk
source activate /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/slamdunk

slamdunk all 	$output_file -r $genome -b $bed_file -o $output_directory -t $cores  -rl $read_length  

alleyoop rates $output_directory/filter/$output_file'_slamdunk_mapped_filtered.bam' -o $output_directory -r $genome # -t $cores

alleyoop utrrates $output_directory/filter/$output_file'_slamdunk_mapped_filtered.bam' -o $output_directory -r $genome  -b $bed_file -l $read_length #-t $cores

source deactivate
##########################


############# Additional plots
source activate r-343
R --slave     --args  $output_directory/$output_file'_slamdunk_mapped_filtered_mutationrates_utr.csv'  $output_directory $output_file /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343 $max_ylim < /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/SLAM-dunk_Global_rate_plotter.R

R --slave     --args  $output_directory/$output_file'_slamdunk_mapped_filtered_overallrates.csv'  $output_directory $output_file /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343 $max_ylim < /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/SLAM-dunk_Overall_rate_plotter.R

source deactivate
##############################


############### Track
source activate r-343

R       --slave --args  $output_directory        $output_file'_complete'     $output_directory'/map/'$output_file'_slamdunk_mapped.bam'        1       200     $read_length     FALSE   transcriptomicsStranded $cores  none    /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343  10      <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/smoothRle.R

R       --slave --args  $output_directory        $output_file'_filtered'     $output_directory'/filter/'$output_file'_slamdunk_mapped_filtered.bam'        1       200     $read_length     FALSE   transcriptomicsStranded $cores  none    /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343  10      <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/smoothRle.R

source deactivate
####################



rm $output_file
rm Rplots.pdf
rm *_fastqc.zip

exit
