commandArgs(FALSE)

peaks =         commandArgs()[4]
directory =	commandArgs()[5]                # home directory of the ChIP not the one of the peak calling
tag =           commandArgs()[6]                # prefix of the file generated
percIdentity =  commandArgs()[7]                # [80 | 85 | 90] for mm and  hs, standard = 85
matrices =      commandArgs()[8]                # [Jaspar | Selex]
minMatches =    as.numeric(commandArgs()[9])    # this indicates to remove the motifs that have less than #matches on the peaks genome wide, standard = 80
foldTreshold =	as.numeric(commandArgs()[10])	# standard = 5
cores =         as.numeric(commandArgs()[11])
lib =           commandArgs()[12]



.libPaths(c(.libPaths(), lib))
library(motifChIPseq)

identity =      paste(percIdentity, "%", sep="")
genome <- Hsapiens
genChar <- "Hsapiens"


############################################################################
param <- as.list(.GlobalEnv)

param$directory <- file.path(param$directory, paste(param$tag, param$matrices, param$percIdentity, sep="_"))
if (! file.exists(param$directory)) dir.create(param$directory)
setwd(param$directory)

pwmList <- get(load(paste(param$lib, "/motifChIPseq/data/", tolower(param$matrices), "Pwm.rda", sep="")))
rm(list = setdiff(ls(), c("param", "pwmList")))


peaks <- import.bed(param$peaks)

#
motScan <- unlist(GRangesList(
mclapply(seq_along(pwmList), function(i) motifScan(pwm = pwmList[[i]], ranges = peaks, genome = param$genome, identity = param$identity, name = names(pwmList)[i]), mc.cores = param$cores)
))

#
motScan <- unlist(GRangesList(  mclapply(chunk(seq_along(motScan), param$cores * 10), function(i) overlapsPeaksMotifs(peaks, motScan[i]), mc.cores = param$cores)  ))
mcols(motScan)$corrDist <- (start(motScan) + width(motScan)/2 - mcols(motScan)$peak_start) / (mcols(motScan)$peak_end  - mcols(motScan)$peak_start)

#
summary_table <- table(mcols(motScan)$peak_name, mcols(motScan)$name)

summary_table_binary <- as.numeric(summary_table > 0)
dim(summary_table_binary) <- dim(summary_table)
colnames(summary_table_binary) <- colnames(summary_table)



#
da <- data.frame(motif_matches = NA, sum_genomic =  get(load(paste(param$lib, "/motifChIPseq/data/", "sumGenomic", param$matrices, "_",  param$genChar, "Motifs", param$percIdentity, ".rda", sep=""))))

da$motif_matches <- as.numeric(table(factor(mcols(motScan)$name, levels = rownames(da))))
da <- motifEnrichment(da$motif_matches, da$sum_genomic, motif_name = rownames(da))

ss <- which(p.adjust(da$p_nbinomial, method="BH") < 0.05 & da$enrichment_ratio > param$foldTreshold)


jpeg("Dataset_vs_total_genomic.jpg", width = 800, height = 800)
  par(mar = c(5,5,4,4))
  plot(da$sum_genomic + 1, da$motif_matches + 1, log = "xy", ylab = "dataset matches", xlab = "genomic matches (total)", pch = 19, cex = 1, cex.lab = 2)
dev.off()


dir.create("Relative_position_within_peaks")
lapply (da[da$enrichment_ratio > param$foldTreshold, "motif_name"], function(nm) {
    if(sum(mcols(motScan)$name == nm) > 2) {
        jpeg(paste0("Relative_position_within_peaks/", nm, ".jpg"), width = 400, height = 400)
        try(plot(density(mcols(motScan)$corrDist[mcols(motScan)$name == nm]), main = paste(nm, "\n", sum(mcols(motScan)$name == nm), " matches"), col = "blue", lwd = 2))
        dev.off()
    }
})

jpeg("Relative_position_within_peaks/All_motifs.jpg", width = 400, height = 400)
	try(plot(density(mcols(motScan)$corrDist), col = "grey", lwd = 2, main = "All motifs"))
dev.off()


save(da, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_da.RData"))
save(motScan, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_motScan.RData"))
save(summary_table, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_summary_table.RData"))
save.image(paste(param$directory, "/Motifs", param$tag, param$matrices, param$percIdentity, ".RData", sep=""))
