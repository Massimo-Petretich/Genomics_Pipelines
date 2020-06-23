commandArgs(FALSE)

directory <-	commandArgs()[4]
tag <-		commandArgs()[5]
bamFile <-	commandArgs()[6]
window <-	as.numeric(commandArgs()[7])	#size of the sliding window
resolution <-	as.numeric(commandArgs()[8])	#bp of the bins used to condense the tracks
fragmentSize <-	as.numeric(commandArgs()[9])	#the bp to add to the reads to complete the theoretical DNA fragment
smooth <-	as.logical(commandArgs()[10])	#logical. Whether to smooth the track or not.
trackType <-	commandArgs()[11]               #types= "input" | "IP" | "ratio" | "transcriptomicsStranded" | "transcriptomicsUnstranded"
cores <-	as.numeric(commandArgs()[12])
bamFileIN <-	commandArgs()[13]
lib <-		commandArgs()[14]
logBase <-	as.numeric(commandArgs()[15])


.libPaths(c(.libPaths(), lib))

library(MPgraphics)
library(parallel)
library(GenomicAlignments)


######################################################################
if (trackType == "IP" | trackType == "input" ) {

	cvg <- MPgraphics:::extractCoverageFromBAM(file = bamFile, fragmentSize = fragmentSize)

	system.time(
		x <- parallel:::mclapply(
			cvg[ sapply(cvg, length) > 3 * ceiling(window/(2*resolution))*resolution ],
			MPgraphics:::smoothRle, 
			window = window, 
			resolution = resolution, 
			sumAggregate = FALSE, 
			smooth = smooth, 
			mc.cores = cores
		)
	)
	x <- IRanges:::RleList(x, compress=FALSE)
	rtracklayer:::export.bw(x, paste(directory, "/", tag,".bw", sep = ""))
}
#########################################################################




#########################################################################
if (trackType == "ratio")   	{
	
	cvg <- MPgraphics:::extractCoverageFromBAM(file = bamFile, fragmentSize = fragmentSize)
	
	system.time(
		x <- parallel:::mclapply(
			cvg[ sapply(cvg, length) > 3 * ceiling(window/(2*resolution))*resolution ], 
			MPgraphics:::smoothRle, 
			window = window, 
			resolution = resolution, 
			sumAggregate = FALSE, 
			smooth = smooth, 
			mc.cores = cores
		)
	)
	
	cvgin <- MPgraphics:::extractCoverageFromBAM(file = bamFileIN, fragmentSize = fragmentSize)

	system.time(
		y <- parallel:::mclapply(
			cvgin[ sapply(cvgin, length) > 3 * ceiling(window/(2*resolution))*resolution ], 
			MPgraphics:::smoothRle, 
			window = window, 
			resolution = resolution, 
			sumAggregate = FALSE, 
			smooth = smooth, 
			mc.cores = cores
		)
	)

	y <- y[names(y) %in% names(x)]

	x <- x[names(x) %in% names(y)]

	system.time(
		z <- parallel:::mcmapply(
			MPgraphics:::logRatio,
			x,
			y, 
			standardise = "mad", 
			mc.cores = cores
		)
	)
	
	z <- IRanges:::RleList(z, compress = FALSE)
	
	rtracklayer:::export.bw(z, 
		paste(directory, "/", tag,"_ratio.bw", sep="")
	)
}
################################################################################




#################################################################################
if (trackType=="both")    {
	
	cvg <- MPgraphics:::extractCoverageFromBAM(file=bamFile, fragmentSize=fragmentSize)
    	
	system.time(
		x <- parallel:::mclapply(
			cvg[ sapply(cvg, length) > 3 * ceiling(window/(2*resolution))*resolution ], 
			MPgraphics:::smoothRle, 
			window = window, 
			resolution = resolution, 
			sumAggregate = FALSE, 
			smooth = smooth, 
			mc.cores = cores
		)
	)
	
	cvgin <- MPgraphics:::extractCoverageFromBAM(file = bamFileIN, fragmentSize = fragmentSize)

	system.time(
		y <- parallel:::mclapply(
			cvgin[ sapply(cvgin, length) > 3 * ceiling(window/(2*resolution))*resolution ], 
			MPgraphics:::smoothRle, 
			window = window, 
			resolution = resolution,  
			sumAggregate = FALSE, 
			smooth = smooth, 
			mc.cores = cores
		)
	)
	
	y <- y[names(y) %in% names(x)]
	
	x <- x[names(x) %in% names(y)]
	
	system.time(
		z <- parallel:::mcmapply(
			MPgraphics:::logRatio,
			x,
			y, 
			standardise = "mad", 
			mc.cores = cores
		)
	)
	
	z <- IRanges:::RleList(z, compress = FALSE)

	rtracklayer:::export.bw(z, paste(directory, "/", tag,"_ratio.bw", sep=""))

	x <- IRanges:::RleList(x, compress = FALSE)

	rtracklayer:::export.bw(x, paste(directory, "/", tag,".bw", sep=""))

}
#########################################################################


#########################################################################
if (trackType == "transcriptomicsStranded" | trackType == "transcriptomicsUnstranded") {


extractCoverageFromBAM2 <- function (file, ...) {

        argnames <- names(list(...))
        addArgs <- list(...)

        if ("range" %in% argnames)      aln1 <- readGAlignmentsFromBam(file, param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = F), which = addArgs$range))
        else                            aln1PlusMinus <- readGAlignments(file,  param=ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = F)))

        aln1Plus <- aln1PlusMinus[strand(aln1PlusMinus)=="+"]
        aln1Minus <- aln1PlusMinus[strand(aln1PlusMinus)=="-"]

        readSplitterMinus <- function(aln1) {
                aln1Whole <- aln1[! grepl(pattern = "N", x = cigar(aln1))]
                aln1Junctions <- aln1[grep(pattern = "N", x = cigar(aln1))]
                aln1Cigar <- cigar(aln1Junctions)
                aln1Whole <- as(aln1Whole, "GRanges")
                aln1Junctions <- as(aln1Junctions, "GRanges")
                aln1CigarsSizes <- lapply(aln1Cigar, function(x) as.integer(strsplit(x, split = "[A-Z]")[[1]]))
                aln1CigarsLetters <- lapply(aln1Cigar, function(x) {try <- strsplit(x, split = "[0-9]")[[1]]; try[try!=""]})

                aln1Junctions <- aln1Junctions[sapply(aln1CigarsLetters, function(x) length(x)==3)]
                aln1CigarsSizes <- aln1CigarsSizes[sapply(aln1CigarsLetters, function(x) length(x)==3)]
                aln1CigarsLetters <- aln1CigarsLetters[sapply(aln1CigarsLetters, function(x) length(x)==3)]

                aln1Junctions <- aln1Junctions[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsSizes <- aln1CigarsSizes[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsLetters <- aln1CigarsLetters[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsSizes <- do.call("rbind", aln1CigarsSizes)

                aln1Junctions2 <- aln1Junctions
                end(aln1Junctions) <- start(aln1Junctions) + aln1CigarsSizes[,1] - 1
                start(aln1Junctions2) <- end(aln1Junctions2) - aln1CigarsSizes[, 3] + 1
                aln2 <- c(aln1Whole, aln1Junctions, aln1Junctions2)
        }

        aln2Minus <- readSplitterMinus(aln1Minus)

        readSplitterPlus <- function(aln1) {
                aln1Whole <- aln1[! grepl(pattern = "N", x = cigar(aln1))]
                aln1Junctions <- aln1[grep(pattern = "N", x = cigar(aln1))]
                aln1Cigar <- cigar(aln1Junctions)
                aln1Whole <- as(aln1Whole, "GRanges")
                aln1Junctions <- as(aln1Junctions, "GRanges")
                aln1CigarsSizes <- lapply(aln1Cigar, function(x) as.integer(strsplit(x, split = "[A-Z]")[[1]]))
                aln1CigarsLetters <- lapply(aln1Cigar, function(x) {try <- strsplit(x, split = "[0-9]")[[1]]; try[try!=""]})

                aln1Junctions <- aln1Junctions[sapply(aln1CigarsLetters, function(x) length(x)==3)]
                aln1CigarsSizes <- aln1CigarsSizes[sapply(aln1CigarsLetters, function(x) length(x)==3)]
                aln1CigarsLetters <- aln1CigarsLetters[sapply(aln1CigarsLetters, function(x) length(x)==3)]

                aln1Junctions <- aln1Junctions[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsSizes <- aln1CigarsSizes[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsLetters <- aln1CigarsLetters[sapply(aln1CigarsLetters, function(x) all(x == c("M","N","M")))]
                aln1CigarsSizes <- do.call("rbind", aln1CigarsSizes)

                aln1Junctions2 <- aln1Junctions
                end(aln1Junctions) <- start(aln1Junctions) + aln1CigarsSizes[,1] -1
                start(aln1Junctions2) <- end(aln1Junctions2) - aln1CigarsSizes[, 3] + 1
                aln2 <- c(aln1Whole, aln1Junctions, aln1Junctions2)
        }

        aln2Plus <- readSplitterPlus(aln1Plus)

        if("type" %in% argnames) {
                if(addArgs$type == "transcriptomicsStranded") {
                        cvgPlus <- coverage(aln2Plus)
                        cvgMinus <- coverage(aln2Minus)
                        cvg <- list(cvgPlus=cvgPlus, cvgMinus=cvgMinus)
                }
                else {
                        aln2 <- c(aln2Plus, aln2Minus)
                        cvg <- coverage(aln2)
                }
        }
        else {
                aln2 <- c(aln2Plus, aln2Minus)
                cvg <- coverage(aln2)
        }

        return(cvg)
}


	if(trackType == "transcriptomicsStranded") {
		cvg <- extractCoverageFromBAM2(file = bamFile,  type = trackType)
		rtracklayer:::export.bw(cvg[[2]], paste(directory, "/", tag,"_Minus.bw", sep = ""))
		rtracklayer:::export.bw(cvg[[1]], paste(directory, "/", tag,"_Plus.bw", sep = ""))
	}
	else {
		cvg <- extractCoverageFromBAM2(file = bamFile, type = trackType)
		rtracklayer:::export.bw(cvg, paste(directory, "/", tag,"_Unstranded.bw", sep = ""))
	}	

}
#########################################################################






# save.image(paste(directory, "/", tag,".RData", sep="") ) 
