commandArgs(FALSE)

queryRegionsFileGR <-	commandArgs()[4]	# Path for a GRanges object saved as .RData	
bamFilesGR <-           commandArgs()[5]	# .txt file containing the paths of .bam files        
fileName <-             commandArgs()[6]
cores <-                as.numeric(commandArgs()[7])



library(rtracklayer)
library(GenomicAlignments)
library(parallel)


loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

queryRegionsGR <- loadRData(queryRegionsFileGR)

bamFilesGR <- readLines(bamFilesGR)

bamGR <- mclapply(bamFilesGR,  function(x) mget(load(x))[[1]], mc.cores = cores)

names(bamGR) <- sapply(bamFilesGR, function(x) {x = basename(x); gsub("_bam_to_GR.RData", "", x)})

tabCount <- mclapply(bamGR, function(x) assay(summarizeOverlaps(queryRegionsGR, x)), mc.cores = cores)
tabCount <- do.call("cbind", tabCount)
colnames(tabCount) <- names(bamGR)

fpkm <- do.call("cbind", lapply(colnames(tabCount),  function(col) tabCount[, col] / (length(bamGR[ names(bamGR) == col ] [[1]]) / 1e6) / (width(queryRegionsGR) / 1e3)) )
colnames(fpkm) <- paste0(colnames(tabCount), "_fpkm")
tabCount <- cbind(tabCount, fpkm)

mcols(queryRegionsGR) <- cbind(mcols(queryRegionsGR), as.data.frame(tabCount))

save(queryRegionsGR, file = fileName)
# save.image("/GWD/bioinfo/projects/RD-Cellzome_NGS/Analysis_Massimo_Petretich/tests/ChIP_table_ws.RData")
