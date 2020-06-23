commandArgs(FALSE)

gtfFile <-		commandArgs()[4]	# path to the .gtf file to use as annotation
bamFiles <-		commandArgs()[5]        # .txt file containing one .bam file path per line
fileName <-             commandArgs()[6]	# .RData file
cores <-                as.numeric(commandArgs()[7])
strandSpecific <-	as.integer(commandArgs()[8])	# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
GTFattrType <-          commandArgs()[9]
lib <- 			commandArgs()[10]

.libPaths(c(.libPaths(), lib))

library(rtracklayer)
library(parallel)
library(Rsubread)



bamFiles <- readLines(bamFiles)

counts <- featureCounts(files = bamFiles, annot.ext = gtfFile, isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType = GTFattrType, strandSpecific = strandSpecific, nthreads = cores)

save(counts, file = fileName)
