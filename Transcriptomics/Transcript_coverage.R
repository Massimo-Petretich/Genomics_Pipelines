commandArgs(FALSE)

directory <-  commandArgs()[4]
tag <-	      commandArgs()[5]
gtfTranscriptsFile <-    commandArgs()[6]
fileBw <-     commandArgs()[7]
fileBwPlus <- commandArgs()[8]
fileBwMinus <-commandArgs()[9]
cores <-      as.numeric(commandArgs()[10])
type <-       commandArgs()[11]
lib <-        commandArgs()[12]

.libPaths(c(.libPaths(), lib))
library(MPgraphics)

setwd(directory)


##### coverage
if (type == "transcriptomicsUnstranded") {
  rleLst <- import.bw(fileBw, as = "RleList")
}


if (type == "transcriptomicsStranded") {
  rleLstPlus <- import.bw(fileBwPlus, as = "RleList")
  rleLstMinus <- import.bw(fileBwMinus, as = "RleList")
  rleLstPlus <- rleLstPlus[names(rleLstPlus) %in% names(rleLstMinus)]
  rleLstMinus <- rleLstMinus[names(rleLstMinus) %in% names(rleLstPlus)]
  rleLst <- lapply(unique(c(names(rleLstPlus), names(rleLstMinus))), function(nm) rleLstPlus[[nm]] + rleLstMinus[[nm]])
  names(rleLst) <- unique(c(names(rleLstPlus), names(rleLstMinus)))
  rleLst <- RleList(rleLst, compress = F)
  rm(list = c("rleLstPlus", "rleLstMinus"))
}

##############



transcripts <- import.gff(gtfTranscriptsFile)
transcripts2 <- transcripts[as.numeric(mcols(transcripts)$FPKM) > 0]
transcripts3 <- transcripts2[mcols(transcripts2)$type == "exon"]
transcripts3 <- transcripts3[seqnames(transcripts3) %in% names(rleLst)]
transcripts3 <- transcripts3[as.numeric(mcols(transcripts3)$FPKM) > quantile(as.numeric(mcols(transcripts3)$FPKM), 0.01)]
seqlevels(transcripts3) <- names(rleLst)
transcripts4 <- split(transcripts3, mcols(transcripts3)$gene_id)

rm(list = c("transcripts", "transcripts2", "transcripts3"))




transcriptLengths <- sapply(transcripts4, function(x) sum(width(x)))
transcriptExp <- sapply(transcripts4, function(x) mcols(x)$FPKM[1])
lengthsCounts <- data.frame(geneName = names(transcriptLengths), lengths = transcriptLengths, expression = as.numeric(transcriptExp), stringsAsFactors = F)
lengthsCounts$expressionCut <- cut(x = lengthsCounts$expression, breaks = quantile(lengthsCounts$expression, c(0, 0.2, 0.4, 0.6, 0.8, 1)), include.lowest = T)
lengthsCounts$lengthsCut <- cut(x = lengthsCounts$lengths, breaks = quantile(lengthsCounts$lengths, c(0, 0.2, 0.4, 0.6, 0.8, 1)), include.lowest = T)
lengthsCounts$category <- paste0("Exp_", as.numeric(lengthsCounts$expressionCut), "-", "Len_", as.numeric(lengthsCounts$lengthsCut))



gRangesListCut <- lapply(levels(as.factor(lengthsCounts$category)), function(x) transcripts4[lengthsCounts[lengthsCounts$category == x, "geneName"]])
names(gRangesListCut) <- levels(as.factor(lengthsCounts$category))


###########################
condensateTranscripts <- function(gRangesList, rleLst) {
  condensedTranscripts <- lapply(gRangesList, function(gRanges) {
    condensedTranscript <- Reduce("c",
                                  lapply(gRanges, function(gr) {
                                    rle <- rleLst[gr][[1]]
                                    return(rle)
                                  }
                                  )
    )
    if(as.character(strand(gRanges[1])) == "-") condensedTranscript <- rev(condensedTranscript)
    return(condensedTranscript)
  }
  )
  
  singleBaseCondensedTranscripts <- lapply(condensedTranscripts, function(x) smoothRle(x, resolution = 1, vector = T))
  singleBasePercenileCondensedTranscripts <- lapply(singleBaseCondensedTranscripts, function(x) {names(x) <- as.numeric(names(x)) / length(x); return(x)})
  
  binnedSingleBasePercenileCondensedTranscripts <- lapply(singleBasePercenileCondensedTranscripts, function(x) {y <- sapply(split(x, cut(as.numeric(names(x)), quantile(as.numeric(names(x)), seq(from = 0, to = 1, by = 0.01)), include.lowest = T)), mean)
                                                                                                                names(y) <- 1:100
                                                                                                                y <- y / sum(y)
                                                                                                                return(y)
  }
  ) 
  
  binnedSingleBasePercenileCondensedTranscripts <- do.call("rbind", binnedSingleBasePercenileCondensedTranscripts)
  
  return(binnedSingleBasePercenileCondensedTranscripts)
}
###########################

condensedCoverage <- condensateTranscripts(gRangesListCut[[13]], rleLst)
boxplot(condensedCoverage, outline = F, bty = "n", whisklty = 1, whiskcol = "grey80", col = "grey40", medcol = "white", medlwd = 2, boxlty = 1, boxcol = "grey40", names = F, frame = F, ylab = "Coverage proportion", xlab = "Transcript position")
###################

library(parallel)
condensedCoverage <- mclapply(gRangesListCut, condensateTranscripts, rleLst = rleLst, mc.cores = cores)

dir.create(paste0(directory, "/Coverage_plots"))
lapply(names(condensedCoverage), function(nm) {
  jpeg(paste0(directory, "/Coverage_plots/", nm, ".jpg"), width = 400, height = 200, quality = 100)
  par(mar = c(3,3,2,0), mgp = c(2,1,0))
  boxplot(condensedCoverage[[nm]], outline = F, bty = "n", whisklty = 1, whiskcol = "grey80", col = "grey40", medcol = "white", medlwd = 2, boxlty = 1, boxcol = "grey40", names = F, frame = F, main = paste0(nm, "; ", dim(condensedCoverage[[nm]])[1], " transcripts"), ylab = "Coverage proportion", xlab = "Transcript position: 5' --> 3'", col.main="grey35", col.lab="grey35", col.axis = "grey35")
  dev.off()
  })  


jpeg(paste0(directory, "/Coverage_plots/Expression_histogram.jpg"), width = 400, height = 300, quality = 100)
  hist(log(lengthsCounts$expression, 10), br = 100, main = "", xlab = "log10 Transcript expression (rpkm)")
  abline(v = log(c(quantile(lengthsCounts$expression, c(0, 0.2, 0.4, 0.6, 0.8, 1))), 10), col = "red")
dev.off()

jpeg(paste0(directory, "/Coverage_plots/Length_histogram.jpg"), width = 400, height = 300, quality = 100)
  hist(log(lengthsCounts$lengths, 10), br = 100, main = "", xlab = "log10 Transcript length (bp)")
  abline(v = log(c(quantile(lengthsCounts$lengths, c(0, 0.2, 0.4, 0.6, 0.8, 1))), 10), col = "blue")
dev.off()



  
save.image("./Workspace_cufflink_transcripts.RData")
