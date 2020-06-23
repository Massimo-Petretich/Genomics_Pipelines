#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/cellranger_count_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/cellranger_count_err.txt

fastq_path=''
output_directory=''
sample_prefix=''
output_tag=''
if [ "$output_tag" == '' ]; then output_tag="$sample_prefix"; fi

genome='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/refdata-cellranger-hg19-3.0.0'
expect_cell=5000

cores=16
mem_per_cpu=6

if [ ! -d $output_directory ]; then mkdir $output_directory; fi
cd $output_directory

/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/cellranger-3.0.2/cellranger count \
			--id=$output_tag \
                   	--transcriptome=$genome  \
                   	--fastqs=$fastq_path \
                   	--sample=$sample_prefix \
                   	--expect-cells=$expect_cell \
                   	--localcores=$cores \
                   	--localmem=$(($cores * $mem_per_cpu))


exit
