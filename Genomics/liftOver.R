commandArgs(trailingOnly = FALSE)
directory 	<-	commandArgs()[4]
tag 		<- 	commandArgs()[5]
cores 		<- 	commandArgs()[6]
genome 		<-	commandArgs()[7]
lib 		<-	commandArgs()[8]

library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(parallel)


pc<-import(paste(directory, "/", tag, "_peaks.bed", sep=""))

chains <- mclapply(list.files("/g/spitz/public/generalFiles/liftOver", full.names=TRUE)[grep(genome, list.files("/g/spitz/public/generalFiles/liftOver"))], import.chain, mc.cores=cores)
names(chains) <- sub(".over.chain","",sub(paste(genome,"To", sep=""), "",  list.files("/g/spitz/public/generalFiles/liftOver")[grep(genome, list.files("/g/spitz/public/generalFiles/liftOver"))]))

if (length(chains)>=1) {
	if (!file.exists(file.path(directory, "liftOver"))) {dir.create(file.path(directory, "liftOver"))}
	
	lftOvr<-mclapply(	chains, 
						function (chain) {
							x<-unlist(liftOver(pc, chain))
							names(x)<-seq_len(length(x)); return(x)
						},  
						mc.cores=cores)

	mcmapply(	function(peaks, names) { export.bed(peaks, paste(directory, "/liftOver/", tag, "_", names, ".bed", sep=""))}, 
				peaks=lftOvr, 
				names=names(lftOvr), 
				mc.cores=cores)
}
