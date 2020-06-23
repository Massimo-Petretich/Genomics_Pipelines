commandArgs(FALSE)

peaks =         commandArgs()[4]
directory =	commandArgs()[5]                # home directory of the ChIP not the one of the peak calling
tag =           commandArgs()[6]                # prefix of the file generated
percIdentity =  commandArgs()[7]                # [80 | 85 | 90] for mm and  hs, standard = 85
matrices =      commandArgs()[8]                # fileneme in the package data folder without .rda
cores =         as.numeric(commandArgs()[9])
lib =           commandArgs()[10]



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

print(matrices)

pwmList <- get(load(paste(param$lib, "/motifChIPseq/data/", param$matrices, ".rda", sep="")))
#pwmList <- get(data(matrices))
rm(list = setdiff(ls(), c("param", "pwmList")))
print(pwmList)

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
#
motifMatches = table(factor(mcols(motScan)$name, levels = unique(mcols(motScan)$name)))






save(motifMatches, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_motifMatches.RData"))
save(motScan, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_motScan.RData"))
save(summary_table, file = paste0(param$tag, "_", param$matrices, "_", param$percIdentity, "_summary_table.RData"))
