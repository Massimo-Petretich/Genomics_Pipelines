############# ONLY FOR PALINDROMIC SITES ###############
commandArgs(FALSE)


range    	<- commandArgs()[4]
pattern1 	<- commandArgs()[5] 
pattern2 	<- commandArgs()[6] 
BSgenome   	<- commandArgs()[7]

directory 	<- commandArgs()[8]
cores		<- commandArgs()[9]
lib 		<- commandArgs()[10]

.libPaths(c(.libPaths(), lib))

library(fragmentsFinder)
library(BSgenome, character.only=TRUE)

genome	<- eval(parse(text=BSgenome))





############# ########################## ###############

if(! file.exists(directory)) {dir.create(directory)}
directory <- paste0(directory, "/", chartr(":,","-.", range))
if(! file.exists(directory))  {dir.create( directory )}


fragments <- fragmSelector(range = range, pattern1 = pattern1, pattern2 = pattern2, genome = genome)
fragments <- readPrimerFinder(fragments = fragments, pattern1 = pattern1, pattern2 = pattern2, genome = genome, primerLength = 24)

fragmentsTmAdjLst <- 	  lapply(c(52, 55, 58, 61), primerTmAdjust, fragments = fragments, clippingCycles = 6)
fragmentsTmAdjDiffTemp <- do.call("rbind", fragmentsTmAdjLst)
fragmentsTmAdjDiffTemp <- fragmentsTmAdjDiffTemp[order(fragmentsTmAdjDiffTemp$fragmentName), ]
fragmentsTmAdjDiffTemp <- fragmentsTmAdjDiffTemp[! duplicated(fragmentsTmAdjDiffTemp), ]

###############
pdf(paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), "_fragments_lengths.pdf"), height = 5, width = 2)
	par(mar = c(8,4,4,4))
	boxplot(fragmentsTmAdjDiffTemp [! duplicated(fragmentsTmAdjDiffTemp$fragmentName), c("pattern1Fragment", "fragmentLength")], ylim = c(0,1500), las = 2)
	abline(h = 150, col = "red")
	abline(h = 300, col = "green")
dev.off()

capture.output(	print("Pattern1Fragment"), 
		summary(fragmentsTmAdjDiffTemp [! duplicated(fragmentsTmAdjDiffTemp$fragmentName), c("pattern1Fragment", "fragmentLength")]$pattern1Fragment), 
		print("fragmentLength"), 
		summary(fragmentsTmAdjDiffTemp [! duplicated(fragmentsTmAdjDiffTemp$fragmentName), c("pattern1Fragment", "fragmentLength")]$fragmentLength), 
		file = paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), "_fragments_lengths.txt")
)
###############

fragmentsTmAdjDiffTempFilt <- filterFragments(fragmentsTmAdjDiffTemp)



if (nrow(fragmentsTmAdjDiffTempFilt ) > 0) {
  
  fragSum2MsMtchTmAdjDiffTemp <- 	fragmentSummaryPart2(fragmentsTmAdjDiffTempFilt, genome, cores=cores)
  fragSum2MsMtchTmAdjDiffTemp <- 	fragmentSummaryPart2(fragSum2MsMtchTmAdjDiffTemp , genome, max.mismatch=2, cores=cores) 
  fragSum2MsMtchTmAdjDiffTempFilt <- 	filterFragments(fragSum2MsMtchTmAdjDiffTemp)
  fragSum2MsMtchTmAdjDiffTempFilt <-    fragmentSummaryPart2(fragSum2MsMtchTmAdjDiffTemp , genome, max.mismatch=3, cores=cores) 
 
  # trim the primer to the last 10bp and count the perfect matches
  primerPattern1_last_10bp <- data.frame(primerPattern1 = sapply(strsplit(fragSum2MsMtchTmAdjDiffTempFilt$primerPattern1, split = "*"), function(x) paste(x[(length(x) - 9):length(x)], collapse = "")), stringsAsFactors=F)
  rownames(primerPattern1_last_10bp) <- rownames(fragSum2MsMtchTmAdjDiffTempFilt)
  primerPattern1_last_10bp <- cbind(fragSum2MsMtchTmAdjDiffTempFilt[, 1:6], primerPattern1_last_10bp)
  try <- fragmentSummaryPart2(primerPattern1_last_10bp, genome, cores=cores)
  fragSum2MsMtchTmAdjDiffTempFilt$genomicMatches_10bp_3_prime <- try$genomicMatches
  rm(try)
	
  fragSum2MsMtchTmAdjDiffTempFilt <- getPCRtemplate(fragSum2MsMtchTmAdjDiffTempFilt, pattern1, range, genome, flankingPattern1=5)
  
  # write the complete csv file
  write.table(	x = fragSum2MsMtchTmAdjDiffTempFilt, 
		file = paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), ".csv"), 
		quote = F, 
		sep = ",", 
		row.names = FALSE
  )

  # print statistics 
  capture.output( print("Filtering for size:"),
                  fragmentsTmAdjDiffTempFilt <- filterFragments(fragmentsTmAdjDiffTemp),
                  print("Filtering for genomic matches:"),
                  fragSum2MsMtchTmAdjDiffTempFilt <- filterFragments(fragSum2MsMtchTmAdjDiffTemp), 
                  file = paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), "_filtering_summary.txt")
  )
  
  # write bed file  
  write.table(fragSum2MsMtchTmAdjDiffTempFilt[!duplicated(fragSum2MsMtchTmAdjDiffTempFilt$fragmentName), c("chr", "start", "end", "fragmentName")], 
              file = paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), ".bed"), 
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t"
  )
  
  # save.image(file = paste0(directory, "/", gsub(pattern = ":|,", replacement = "-", x = range), ".RData"))

}

range_whole_chromosome <- paste0(strsplit(range, ":")[[1]][1], ":1-",  length(genome[[which(seqnames(genome) %in% strsplit(range, ":")[[1]][1])]]) )
fragments_whole_chromosome <- fragmSelector(range = range_whole_chromosome, pattern1 = pattern1, pattern2 = pattern2, genome = genome)
save(fragments_whole_chromosome, file = "fragments_whole_chromosome.RData")
