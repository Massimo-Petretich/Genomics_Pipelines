commandArgs(FALSE)

bamFile <-	commandArgs()[4]
directory <- 	commandArgs()[5]
tag <-		commandArgs()[6]		# The prefix for the file name 
fragmentSize <-	as.numeric(commandArgs()[7])	# The total fragment size to adjust each read if extend is set to TRUE
extend <-  	commandArgs()[8]		# Logical



library(rtracklayer)
library(GenomicAlignments)


bam <- readGAlignments(bamFile)

bam <- as(bam, "GRanges")

if(extend) bam <- resize(bam, width = fragmentSize, fix = "start", ignore.strand = F)

save(bam, file = paste0(directory, "/", tag, "_bam_to_GR.RData"))
