commandArgs(FALSE)

curTab = commandArgs()[4]
directory = commandArgs()[5]
tag = commandArgs()[6]
lib = commandArgs()[7] 
max_ylim = as.numeric(commandArgs()[8])

.libPaths(c(.libPaths(), lib))

library(ggplot2)

curTab = read.table(curTab, stringsAsFactors=FALSE)

curTab[, c("A", "C", "G", "T")] <- curTab[, c("A", "C", "G", "T")] / rowSums(curTab[, c("A", "C", "G", "T")]) * 100
curTab[, c("a", "c", "g", "t")] <- curTab[, c("a", "c", "g", "t")] / rowSums(curTab[, c("a", "c", "g", "t")]) * 100

printTab = data.frame(rates=c(rep("AT",2),rep("AC",2),rep("AG",2),
                              rep("TA",2),rep("TC",2),rep("TG",2),
                              rep("CA",2),rep("CT",2),rep("CG",2),
                              rep("GA",2),rep("GT",2),rep("GC",2)), strand = rep(c("+","-"),12),
                      rate_percent = c(curTab["A","T"],curTab["A","t"],curTab["A","C"],curTab["A","c"],curTab["A","G"],curTab["A","g"],
                                       curTab["T","A"],curTab["T","a"],curTab["T","C"],curTab["T","c"],curTab["T","G"],curTab["T","g"],
                                       curTab["C","A"],curTab["C","a"],curTab["C","T"],curTab["C","t"],curTab["C","G"],curTab["C","g"],
                                       curTab["G","A"],curTab["G","a"],curTab["G","T"],curTab["G","t"],curTab["G","C"],curTab["G","c"])
)




printTab$y = -0.3
printTab[printTab$strand == "-", ]$y = printTab[printTab$strand == "-", ]$rate_percent + printTab[printTab$strand == "+", ]$rate_percent

curPlot = qplot(x=rates, y=rate_percent, fill=strand,data=printTab) + 
  ylim(0, max_ylim) + 
  geom_bar(stat="identity") + 
  ylab("Rate percent %") + 
  ggtitle(label = tag) +
  xlab("") +
  theme(plot.title = element_text(color = "grey50", size = 22),
        axis.title.y = element_text(color = "grey50", size = 16),
        axis.text.x = element_text(color = "grey50", size = 16, angle = 90),
        axis.text.y = element_text(color = "grey50", size = 16),
        legend.title = element_text(color = "grey50", size = 16),
        legend.text = element_text(color = "grey50", size = 16)
  )

curPlot + ggsave(paste0(directory, "/", tag, "_overall_mutation_rates.jpg"), device = "jpeg", width = 6, height = 4)
