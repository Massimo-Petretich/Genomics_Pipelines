commandArgs(FALSE)

cores		= as.numeric(commandArgs()[4])
logBase		= as.numeric(commandArgs()[5])
treshold	= as.numeric(commandArgs()[6])
tag		= commandArgs()[7]
pval		= as.numeric(commandArgs()[8])
directory 	= commandArgs()[9]

IPbw  		= commandArgs()[10]
INbw  		= commandArgs()[11]
INbam 		= commandArgs()[12]
IPbam 		= commandArgs()[13]

resolution	= as.numeric(commandArgs()[14])
smoothInput	= as.logical(commandArgs()[15])
lib		= commandArgs()[16]

.libPaths(c(.libPaths(), lib))

logRatio <- function(x,y, logBase){    
        x<-window(x,start=1, end=min(length(x),length(y)))    
        y<-window(y,start=1, end=min(length(x),length(y)))    
	
	medianY <- median(y)
	medianX <- median(x)
	sdY <- mad(y)
	sdX <- mad(x)
	x <- x * (medianY/medianX) * ((sdY/medianY)^(1/16)/(sdX/medianX)^(1/16))

        z <- log(x/y,logBase)
        runValue(z)[ runValue(z) == "-Inf" |  runValue(z) == "Inf" | runValue(z) == "NaN"] <- 0
        z <- z-(sum(runLength(z)*runValue(z))/length(z))
        return(z)
} 


library(rtracklayer)

x<-import.bw(IPbw, as="RleList")

y<-import.bw(INbw, as="RleList")

if (! smoothInput) {
	library(MPgraphics, lib.loc  = lib)
	window <- 2000
	y <- mclapply(y[ sapply(y, length) > 3 * ceiling(window/(2*resolution))*resolution ], smoothRle, window = window, resolution = resolution, sumAggregate = F, smooth = T, mc.cores = cores)
	y <- IRanges:::RleList(y, compress=FALSE)

}
 
y <- y[names(y) %in% names(x)]

x <- x[names(x) %in% names(y)]

z <- mcmapply(logRatio,x,y, logBase = logBase, mc.cores = cores)
z <- RleList(z, compress=FALSE)
zGR <- as(z,"GRanges")

qtl <- mclapply(z, quantile, probs = 1-treshold, mc.cores = cores)
qtl <- sum(sapply(z,length)*(unlist(qtl))) / sum(as.numeric(sapply(z,length)))

zQtl <- zGR[mcols(zGR)[,1] >= qtl]
zQtl <- reduce(zQtl)
seqlengths(zQtl) <- rep(NA, length(seqlengths(zQtl)))





library(Rsamtools)
library(GenomicAlignments)


bamGR <- mclapply(c(INbam, IPbam),  function(x) mget(load(x))[[1]], mc.cores = cores)
names(bamGR) <- c("Input", "IP")

count <- mclapply(bamGR, function(x) assay(summarizeOverlaps(zQtl, x)), mc.cores = cores)
count <- do.call("cbind", count)
colnames(count) <- names(bamGR)


#files <- BamFileList(c(INbam, IPbam))
#count <- assay(summarizeOverlaps(zQtl, files))

fr <- median(x[[1]]/(x[[1]]+y[[1]]), na.rm=TRUE)


p  <- ppois(count[, 2], rowSums(count)*fr, lower.tail = F)
p[rowSums(count) == 0] <- NA
pa <- p.adjust(p, method = "BH")

mcols(zQtl)["score"] <- round(count[, 2] / count[, 1], 2)
mcols(zQtl)["padj"] <- pa
zQtlPA <- zQtl[pa <= pval & !is.na(pa)]
mcols(zQtlPA)["name"] <- paste("LG_peaks_", order(mcols(zQtlPA)[, 1], decreasing = T), sep = "")

export.bed(zQtlPA, paste0(directory, "/", tag, "_peaksLG.bed"))
save(zQtlPA, file = paste0(directory, "/", tag, "_peaksLG.RData"))

#save.image(paste0(directory, "/", tag, "_peaksLG.RData"))


