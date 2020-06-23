#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/RNA-seq_STAR_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/RNA-seq_STAR_err.txt




fastqFile=''
outputDirectory=''
outputFile=''

cores='8'

genome='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/'
annotation='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Annotations/hg19_Ensembl_GRCh37.gtf'

type='transcriptomicsUnstranded'

fastq5primeTrim=13      # The first base to keep from the fastq read
fastq3primeTrim=0       # The nucleotide to remove from the 3 prime end, default 0


fastqMinQuality=30      # Trim reads from quality below this treshold
fastqMinLength=34       # Discard shorter reads

discardFlag=260 	# discard unmapped and secondary alignments from sam file
discardMappingQuality=5 #Discard lower quality than 5 from sam file



################################################################
source activate r-343

if [ ! -d $outputDirectory ];			then  	mkdir $outputDirectory;			fi  
if [ ! -d $outputDirectory/$outputFile ];	then	mkdir $outputDirectory/$outputFile;	fi	

outputDirectory=$outputDirectory/$outputFile
cd $outputDirectory



readLength=$(gunzip -c $fastqFile | head -n 4 | awk 'FNR==2{print $1}' | wc | awk '{print $3}')
totalTrim=$(($fastq5primeTrim + $fastq3primeTrim))
usefulRead=$(($readLength - $totalTrim))

readMismatches=$(($usefulRead / 19))

lastBaseKept=$(expr $readLength - $fastq3primeTrim)



date	>  $outputDirectory/summary.txt
declare -p | awk '{print $3}' | egrep "fastqFile|outputFile|outputDirectory|cores|genome|type|annotation|fastq5primeTrim|fastq3primeTrim|fastqMinQuality|fastqMinLength|discardFlag|discardMappingQuality|readMismatches" >> $outputDirectory/summary.txt



############# ALIGNING #######################
gunzip -c  $fastqFile | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastx_trimmer -f $fastq5primeTrim -l $lastBaseKept -Q33 | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastq_quality_trimmer -t $fastqMinQuality -Q33 -l $fastqMinLength | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastp --stdin  --stdout --trim_poly_g  --trim_poly_x | perl /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/prinseq-lite-0.20.4/prinseq-lite.pl -trim_tail_right 10 -min_len $fastqMinLength -fastq stdin  -out_good $outputDirectory/$outputFile


/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/STAR-2.7.0f/bin/Linux_x86_64/STAR --genomeDir $genome --runThreadN $cores --readFilesIn $outputDirectory/$outputFile.fastq --outFileNamePrefix $outputDirectory/$outputFile --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMismatchNoverLmax 0.06 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 

cat $outputDirectory/$outputFile'Log.final.out' >> $outputDirectory/summary.txt

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools view -F $discardFlag  -q $discardMappingQuality       $outputDirectory/$outputFile'Aligned.sortedByCoord.out.bam' -b > $outputDirectory/$outputFile.bam

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools index $outputDirectory/$outputFile.bam


echo "Original fastq length: "$(( $(gunzip -c $fastqFile  | wc -l)/4)) >> $outputDirectory/summary.txt
echo "Trimmed fastq length: "$(($(cat $outputDirectory/$outputFile.fastq | wc -l)/4)) >> $outputDirectory/summary.txt

echo "Original bam length: "$(samtools view $outputDirectory/$outputFile'Aligned.sortedByCoord.out.bam' | wc -l) >> $outputDirectory/summary.txt
echo "Quality filtered bam length: "$(samtools view $outputDirectory/$outputFile.bam | wc -l) >> $outputDirectory/summary.txt

rm -r $(find $outputDirectory -name "*STARtmp")
################################################


############ FASTQC ##############################################
/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $fastqFile --outdir=$outputDirectory
mv $(basename "$fastqFile" | sed 's/\.txt\.gz/_fastqc\.html/') $outputDirectory/$outputFile'_original_fastqc.html'

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $outputDirectory/$outputFile.fastq   --outdir=$outputDirectory
mv $(echo $outputDirectory/$outputFile'_fastqc.html') $outputDirectory/$outputFile'_trimmed_fastqc.html'
#################################################################



#################### TRACKS ####################
R       --slave         --args          $outputDirectory        $outputFile             $outputDirectory/$outputFile.bam        1       200     $usefulRead  FALSE      $type   $cores  none    /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343    10      <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/smoothRle.R
################################################


for i in $( ls $outputDirectory | egrep -vi "^$outputFile.bam|summary.txt|$outputFile.*.bw|fastqc.html$|.pdf$|.RData$|deletions.bed|junctions.bed|insertions.bed|transcripts.gtf|transcripts.RData|transcripts.bed$" ); do if [ ! -d $outputDirectory/$i ]; then rm -rf $outputDirectory/$i; fi; done

date >>  $outputDirectory/summary.txt



exit
