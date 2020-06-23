#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/RNA-seq_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/RNA-seq_err.txt




fastqFile=''
outputDirectory=''
outputFile=''

cores='8'

genome='/GWD/bioinfo/projects/RD-TSci-CommonData/all_datasets/genomes/genome_hg19/hg19/bowtie2.0.1_index/hg19_complete'
annotation='/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Annotations/hg19_Ensembl_GRCh37.gtf'

type='transcriptomicsUnstranded'

fastq5primeTrim=13      # The first base to keep from the fastq read
fastq3primeTrim=0       # The nucleotide to remove from the 3 prime end, default 0


fastqMinQuality=30      # Trim reads from quality below this treshold
fastqMinLength=34       # Discard shorter reads

discardFlag=260 	# discard unmapped and secondary alignments from sam file
discardMappingQuality=5 #Discard lower quality than 5 from sam file

# segmentLength=18        # For Tophat, length of the splitted read for exon-exon spanning reads (should be ~ 1/2 read length)


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

if [ $usefulRead -gt 50 ] 
        then
        segmentLength=25  
	segmentMismatches=2     
	else
        segmentLength=$(($usefulRead / 2))
	segmentMismatches=$(($segmentLength / 12))
fi

if [ $readMismatches -gt 2 ] 
	then
	readEditDist=$readMismatches
        else
	readEditDist=2
fi

lastBaseKept=$(expr $readLength - $fastq3primeTrim)



date	>  $outputDirectory/summary.txt
declare -p | awk '{print $3}' | egrep "fastqFile|outputFile|outputDirectory|cores|genome|type|annotation|fastq5primeTrim|fastq3primeTrim|fastqMinQuality|fastqMinLength|discardFlag|discardMappingQuality|segmentLength|segmentMismatches|readMismatches" >> $outputDirectory/summary.txt



############# ALIGNING AND QUALITY CONTROL #######################
gunzip -c  $fastqFile | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastx_trimmer -f $fastq5primeTrim -l $lastBaseKept -Q33 | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastq_quality_trimmer -t $fastqMinQuality -Q33 -l $fastqMinLength | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastp --stdin  --stdout --trim_poly_g  --trim_poly_x | perl /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/prinseq-lite-0.20.4/prinseq-lite.pl -trim_tail_right 10 -min_len $fastqMinLength -fastq stdin  -out_good $outputDirectory/$outputFile


# gunzip -c  $fastqFile | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastx_trimmer -f $fastq5primeTrim -Q33 | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastq_quality_trimmer -t $fastqMinQuality -Q33 -l $fastqMinLength -o $outputDirectory/$outputFile.fastq


/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/tophat -p $cores	--segment-mismatches=$segmentMismatches  --segment-length=$segmentLength	--read-mismatches=$readMismatches	--read-edit-dist=$readEditDist	-o $outputDirectory -G $annotation $genome $outputDirectory/$outputFile.fastq   2>> $outputDirectory/summary.txt

mv $outputDirectory/accepted_hits.bam $outputDirectory/NF_$outputFile.bam
/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools view -F $discardFlag  -q $discardMappingQuality       $outputDirectory/NF_$outputFile.bam -b > $outputDirectory/$outputFile.bam

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools index $outputDirectory/$outputFile.bam
##################################################################


echo "Original fastq length: "$(( $(gunzip -c $fastqFile  | wc -l)/4)) >> $outputDirectory/summary.txt
echo "Trimmed fastq length: "$(($(cat $outputDirectory/$outputFile.fastq | wc -l)/4)) >> $outputDirectory/summary.txt

echo "Original bam length: "$(samtools view $outputDirectory/NF_$outputFile.bam | wc -l) >> $outputDirectory/summary.txt
echo "Quality filtered bam length: "$(samtools view $outputDirectory/$outputFile.bam | wc -l) >> $outputDirectory/summary.txt

############ FASTQC ##############################################
/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $fastqFile --outdir=$outputDirectory
mv $(basename "$fastqFile" | sed 's/\.txt\.gz/_fastqc\.html/') $outputDirectory/$outputFile'_original_fastqc.html'

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $outputDirectory/$outputFile.fastq   --outdir=$outputDirectory
mv $(echo $outputDirectory/$outputFile'_fastqc.html') $outputDirectory/$outputFile'_trimmed_fastqc.html'
#################################################################


##############################################################
/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/cufflinks-2.2.1.Linux_x86_64/cufflinks -p  $cores  -o $outputDirectory -g $annotation  $outputDirectory/$outputFile.bam  > /dev/null 2>&1

R       --slave         --args          $outputDirectory/transcripts.gtf	$outputDirectory/$outputFile'_transcripts.bed'	<	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/gtf_to_bed.R

mv $outputDirectory/transcripts.RData $outputDirectory/$outputFile'_transcripts.RData'
##############################################################



#################### TRACKS ####################
readLength=$(/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools view $outputDirectory/$outputFile.bam | head -n 1 | awk '{print $10}' | wc | awk '{print $3}')

R	--slave 	--args  	$outputDirectory	$outputFile		$outputDirectory/$outputFile.bam	1	200	$readLength  FALSE	$type	$cores	none	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343	10	<	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/smoothRle.R
################################################


######### TRANSCRIPT COVERAGE ##################
R       --slave         --args          $outputDirectory        $outputFile             $outputDirectory/transcripts.gtf        $outputDirectory/$outputFile'_Unstranded.bw'    $outputDirectory/$outputFile'_Plus.bw'  $outputDirectory/$outputFile'_Minus.bw' $cores  $type   /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343    <       /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/Transcript_coverage.R
################################################



for i in $( ls $outputDirectory | egrep -vi "^$outputFile.bam|summary.txt|$outputFile.*.bw|.html$|.pdf$|.RData$|deletions.bed|junctions.bed|insertions.bed|transcripts.gtf|transcripts.RData|transcripts.bed$" ); do if [ ! -d $outputDirectory/$i ]; then rm $outputDirectory/$i; fi; done

date >>  $outputDirectory/summary.txt


exit
