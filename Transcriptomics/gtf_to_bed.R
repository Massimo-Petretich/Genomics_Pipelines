commandArgs(FALSE)

gtf <-			commandArgs()[4]        
outputFile <-		commandArgs()[5]



library(rtracklayer)


x <- import.gff(gtf)

mcols(x)$FPKM <- as.numeric(mcols(x)$FPKM)
x2 <- split(mcols(x)$FPKM, mcols(x)$gene_id)
x3 <- sapply(x2, sum)

mcols(x) <- mcols(x)[c("type", "gene_id", "transcript_id", "FPKM")]
x4 <- split(x, mcols(x)$gene_id)
x4 <- x4[x3 > 0]

save(x4, file = gsub("\\.gtf", "\\.RData", gtf))

export.bed(unlist(x4), outputFile)


