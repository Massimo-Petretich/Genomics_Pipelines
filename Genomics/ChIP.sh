#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G 
#SBATCH --time=99:00:00
#SBATCH -o /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/ChIP_out.txt
#SBATCH -e /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tmp/ChIP_err.txt


fastqFile=''

outputFile=''
outputDirectory=''

tagIN='none'
inputDirectory='none'

cores='4'

EXTEND=150
genome='/GWD/bioinfo/projects/RD-TSci-CommonData/all_datasets/genomes/genome_hg19/hg19/bowtie2.0.1_index/hg19_complete'
discardFlag=260	# discard unmapped and secondary alignments from sam file
discardMappingQuality=5	#Discard lower quality than 5 from sam file

type='IP'	# input | IP | ratio | both
window=500
resolution=250
smooth='FALSE'
smoothInput='TRUE'	# For PK_MP if the Input has to be further smoothened

genomeSize='hs'	# For MACS14. Effective genome size. Shortcuts:'hs' for human (2.7e9) 
peakRange='10,30'

MPtreshold=0.1
MPpval=0.05
logBase=10

paired='FALSE'

fastqMinQuality=30
fastqMinLength=45

source activate r-343



if [ "$inputDirectory" == "none" ]; then inputDirectory=$outputDirectory; fi

if [ "$type" == "ratio" ] && [ ! -d "$inputDirectory/$tagIN" ] || [ "$type" == "both" ] && [ ! -d "$inputDirectory/$tagIN" ]; then echo "Specify a valid input control directory or a valid input tag"; exit 1; fi

if [ "$smooth" != "TRUE" ] && [  "$smooth" != "FALSE" ]; then echo "Smooth: logical 'TRUE' | 'FALSE'";	exit 1; fi

if [ "$type" != "input" ] &&  [ "$type" != "IP" ] && [  "$type" != "ratio" ] && [  "$type" != "both" ]; then echo "Specify type 'input' | 'IP' | 'ratio' | 'both'"; exit 1; fi

if [ "$motifs" == "T" ]; then motifs="TRUE"; fi 


if [ ! -d $outputDirectory ];               	then  	mkdir $outputDirectory;            	fi  
if [ ! -d $outputDirectory/$outputFile ];	then	mkdir $outputDirectory/$outputFile;	fi	

outputDirectory=$outputDirectory/$outputFile
cd $outputDirectory

date	>  $outputDirectory/summary.txt
declare -p | awk '{print $3}' | egrep "paired|inputDirectory|fastqFile|outputFile|outputDirectory|cores|EXTEND|genome|type|window|resolution|directoryIN|tagIN|smooth|peakRange|motifs|logBase|MPtreshold|MPpval|fastqMinQuality|fastqMinLength" >> $outputDirectory/summary.txt

echo $HOSTNAME >>  $outputDirectory/summary.txt


############# ALIGNING AND INDEXING #######################
gunzip -c  $fastqFile | /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/fastx_bin/bin/fastq_quality_trimmer -t $fastqMinQuality -Q33 -l $fastqMinLength -o $outputDirectory/$outputFile.fastq

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/bowtie2 -N 1 -p $cores  -x $genome -U $outputDirectory/$outputFile.fastq	-S $outputDirectory/$outputFile.sam	2>> $outputDirectory/summary.txt


/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools view -bS $outputDirectory/$outputFile.sam  > $outputDirectory/NS_$outputFile.bam

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools view -F $discardFlag  -q $discardMappingQuality       $outputDirectory/NS_$outputFile.bam -b > $outputDirectory/F_NS_$outputFile.bam

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools sort -@ $cores -o $outputDirectory/$outputFile.bam $outputDirectory/F_NS_$outputFile.bam

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/samtools index $outputDirectory/$outputFile.bam


echo "Original fastq length: "$(( $(gunzip -c $fastqFile  | wc -l)/4)) >> $outputDirectory/summary.txt
echo "Trimmed fastq length: "$(($(cat $outputDirectory/$outputFile.fastq | wc -l)/4)) >> $outputDirectory/summary.txt

echo "Original bam length: "$(samtools view $outputDirectory/NS_$outputFile.bam | wc -l) >> $outputDirectory/summary.txt
echo "Quality filtered bam length: "$(samtools view $outputDirectory/$outputFile.bam | wc -l) >> $outputDirectory/summary.txt
##########################################################


################# FASTQC #################################
/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $fastqFile                           --outdir=$outputDirectory
mv $(basename "$fastqFile" | sed 's/\.txt\.gz/_fastqc\.html/') $outputDirectory/$outputFile'_original_fastqc.html'

/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/fastqc --java=/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/java $outputDirectory/$outputFile.fastq   --outdir=$outputDirectory
mv $(echo $outputDirectory/$outputFile'_fastqc.html') $outputDirectory/$outputFile'_trimmed_fastqc.html'

# rm -r 	$outputDirectory/*_fastqc/Icons
# rm 		$outputDirectory/*_fastqc/*data.txt
##############################################################



############# PEAK CALLING #################
if [ ! -d $outputDirectory/peakCalling ];   then    mkdir $outputDirectory/peakCalling; fi

if [ $type != input ]
then
	if [   -f $inputDirectory/$tagIN/$tagIN.bam ]
		then
		/GWD/bioinfo/share/scripts/macs14 -t 	$outputDirectory/$outputFile.bam  -c $inputDirectory/$tagIN/$tagIN.bam  	-n $outputDirectory/peakCalling/$outputFile 	-f BAM 	-g $genomeSize -m $peakRange --shiftsize=$(( $EXTEND/2 ))  2>> $outputDirectory/summary.txt
		else
		/GWD/bioinfo/share/scripts/macs14 -t 	$outputDirectory/$outputFile.bam						-n $outputDirectory/peakCalling/$outputFile		-f BAM 	-g $genomeSize -m $peakRange --shiftsize=$(( $EXTEND/2 ))  2>> $outputDirectory/summary.txt
	fi
	
	R --slave --args < $outputDirectory/peakCalling/$outputFile\_model.r
	rm $outputDirectory/peakCalling/*model.r

else
/GWD/bioinfo/share/scripts/macs14 -t    $outputDirectory/$outputFile.bam	-n $outputDirectory/peakCalling/$outputFile -f BAM  -g $genomeSize -m $peakRange --shiftsize=$(( $EXTEND/2 ))  2>> $outputDirectory/summary.txt

R --slave --args < $outputDirectory/peakCalling/$outputFile\_model.r
rm $outputDirectory/peakCalling/*model.r

fi
################################################


#################### TRACKS ####################
if [ "$type" == "ratio" ] && [ -f $inputDirectory/$tagIN/$tagIN.bam ] || [ "$type" == "both" ] && [ -f $inputDirectory/$tagIN/$tagIN.bam ] ||  [ "$type" == "input" ] || [ "$type" == "IP" ]
	then
	R	--slave	--args	$outputDirectory	$outputFile	$outputDirectory/$outputFile.bam	$window	$resolution	$EXTEND	$smooth	$type	$cores	$inputDirectory/$tagIN/$tagIN.bam	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343	$logBase	< /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/smoothRle.R
	else
	echo "INPUT CHROMATIN FILE: $inputDirectory/$tagIN/$tagIN.bam 		Specify a valid Input Chromatin .bam file"
fi
################################################


############### bam to GRanges ################
R     --slave     --args     $outputDirectory/$outputFile.bam   $outputDirectory        $outputFile  $EXTEND     TRUE  < /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/bam_to_GR.R
###############################################



################ PEAK CALLING MP ###############
if [ -f $inputDirectory/$tagIN/$tagIN.bw ] && [ -f $outputDirectory/$outputFile.bw ] && [ -f $inputDirectory/$tagIN/$tagIN'_bam_to_GR.RData' ] && [ -f $outputDirectory/$outputFile'_bam_to_GR.RData' ]
    then
    R     --slave     --args    $cores $logBase $MPtreshold $outputFile   $MPpval $outputDirectory/peakCalling   $outputDirectory/$(ls $outputDirectory | grep .bw | grep -v ratio)   $inputDirectory/$tagIN/$tagIN.bw    $inputDirectory/$tagIN/$tagIN'_bam_to_GR.RData'  $outputDirectory/$outputFile'_bam_to_GR.RData'     $resolution	$smoothInput	/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343	< /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/PK_MP.R
fi

if [ $type == input ]
then	
	MPtresholdInput=0.001
        windowSmall=2000
        windowBig=75000
     R --slave     --args    $cores $logBase $MPtresholdInput $windowSmall $windowBig $outputFile   $outputDirectory/peakCalling   $outputDirectory/$(ls $outputDirectory | grep .bw | grep -v ratio)   $resolution  /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/rlib/r-343    < /GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/Software/Scripts/PK_MP_Input.R
fi
###############################################



if [ ! "$type" == "input" ]
then
	for i in $( ls $outputDirectory | egrep -vi "summary.txt|$outputFile.*.bw|_fastqc$|.pdf$|.*peaks.*|.RData$|.html$" ); do if [ ! -d $outputDirectory/$i ]; then rm $outputDirectory/$i; fi; done
else
	for i in $( ls $outputDirectory | egrep -vi "^$outputFile.bam|summary.txt|$outputFile.*.bw|_fastqc$|.pdf$|.*peaks.*|.RData$|.html$" ); do if [ ! -d $outputDirectory/$i ]; then rm $outputDirectory/$i; fi; done
fi


date >>  $outputDirectory/summary.txt


exit
