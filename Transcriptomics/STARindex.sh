#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/STAR_Index_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/STAR_Index_err.txt

/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/STAR-2.7.0f/bin/Linux_x86_64/STAR --genomeDir /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/ --runMode genomeGenerate --genomeFastaFiles /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-2015-07-17-14-31-42/Genes/genes.gtf --sjdbOverhang 71 --runThreadN 8

exit
