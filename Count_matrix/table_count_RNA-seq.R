commandArgs(FALSE)

queryRegionsFileGR <-	commandArgs()[4]        # Path for a GRanges object saved as .RData
bamFiles <-		commandArgs()[5]        # .txt file containing the paths of .bam files
fileName <-             commandArgs()[6]
cores <-                as.numeric(commandArgs()[7])
ignoreStrand <-		as.logical(commandArgs()[8])
singleEnd <-		as.logical(commandArgs()[9])
mode <- 		commandArgs()[10]

library(rtracklayer)
library(GenomicAlignments)
library(parallel)


loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

queryRegionsGR <- loadRData(queryRegionsFileGR)

bamFiles <- readLines(bamFiles)

count <- mclapply(bamFiles, function(bam) assay(summarizeOverlaps(queryRegionsGR,  BamFileList(bam), mode=mode, ignore.strand=ignoreStrand, singleEnd=singleEnd)), mc.cores = cores)
count <- do.call("cbind", count)
colnames(count) <- gsub(pattern = "\\.bam", "", colnames(count))


save(count, file = fileName)
