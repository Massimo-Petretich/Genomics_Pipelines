commandArgs(FALSE)

annotatedGtf <-		commandArgs()[4]	# path to the .gtf file to use as annotation
mergedGtf <-		commandArgs()[5]        # path to the merged (cuffmerge) .gtf file
fileName <-             commandArgs()[6]	# .RData file
cores <-                as.numeric(commandArgs()[7])

library(rtracklayer)
library(org.Hs.eg.db)
orgDb <- org.Hs.eg.db


annotated <- import.gff(annotatedGtf)
merged <- import.gff(mergedGtf)

mcols(merged)$score <- 0

nm <- mcols(merged)$nearest_ref
nm2 <- split(nm, mcols(merged)$gene_id)
nm2 <- lapply(nm2, unique)

mappings <- data.frame(gene = as.character(mcols(annotated)$gene_id), transcript = as.character(mcols(annotated)$transcript_id), stringsAsFactors=F)
mappings <- mappings[! duplicated(mappings$transcript), ]

nmGene <- Reduce("c", 
                 mclapply(split(nm2, rep(1:cores, each = ceiling(length(nm2) / cores))), 
                   function(x) lapply(x, function(x) unique(mappings[mappings$transcript %in% x, "gene"])), 
                   mc.cores = cores)
)

save(nmGene, file = fileName)
