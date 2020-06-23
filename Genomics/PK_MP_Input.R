commandArgs(FALSE)

cores		= as.numeric(commandArgs()[4])
logBase		= as.numeric(commandArgs()[5])
treshold	= as.numeric(commandArgs()[6])
windowSmall	= as.numeric(commandArgs()[7])
windowBig	= as.numeric(commandArgs()[8])
tag		= commandArgs()[9]
directory 	= commandArgs()[10]

INbw  		= commandArgs()[11]

resolution	= as.numeric(commandArgs()[12])
lib		= commandArgs()[13]

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

y<-import.bw(INbw, as="RleList")

library(MPgraphics, lib.loc  = lib)

# window <- 2000
y <- mclapply(y[ sapply(y, length) > 3 * ceiling(windowSmall/(2*resolution))*resolution ], smoothRle, window = windowSmall, resolution = resolution, sumAggregate = F, smooth = T, mc.cores = cores)
y <- IRanges:::RleList(y, compress=FALSE)

# window <- 100000
y_long_smooth <- mclapply(y[ sapply(y, length) > 3 * ceiling(windowBig/(2*resolution))*resolution ], smoothRle, window = windowBig, resolution = resolution, sumAggregate = F, smooth = T, mc.cores = cores)
y_long_smooth <- IRanges:::RleList(y_long_smooth, compress=FALSE)
 

z <- mcmapply(logRatio, y, y_long_smooth, logBase = logBase, mc.cores = cores)
z <- RleList(z, compress=FALSE)
zGR <- as(z,"GRanges")

qtl <- mclapply(z, quantile, probs = 1-treshold, mc.cores = cores)
qtl <- sum(sapply(z,length)*(unlist(qtl))) / sum(as.numeric(sapply(z,length)))

zQtl <- zGR[mcols(zGR)[,1] >= qtl]
zQtl <- reduce(zQtl)
seqlengths(zQtl) <- rep(NA, length(seqlengths(zQtl)))

export.bed(zQtl, paste0(directory, "/", tag, "_peaksLG.bed"))


# save.image(paste0(directory, "/", tag, "_peaksLG.RData"))
