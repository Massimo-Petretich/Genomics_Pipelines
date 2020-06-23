commandArgs(FALSE)

file_wig <- commandArgs()[4]

library(rtracklayer)

y <- readLines(file_wig)
idx <- grep("variableStep", y)

# empty_idx <- idx[which(idx - c(idx[-1], 0) == -1)]
empty_idx <- idx[which(-1*(idx - c(idx[-1], 0)) < 100)]


if(idx[length(idx)] != length(y)) empty_idx <- empty_idx[-length(empty_idx)]

y2 <- y[-empty_idx]
writeLines(text = y2, con = gsub('.wig', '_filtered.wig', file_wig))

x <- import.wig(gsub('.wig', '_filtered.wig', file_wig))

x_list <- split(x, seqnames(x))

seqlengths(x) <- sapply(x_list, function(x) max(end(x)))

x <- x[-grep('random', as.character(seqnames(x)))]


export.bw(x, con = gsub('.wig', '.bw', file_wig))

