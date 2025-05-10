#' ---
#' title: "PGS - Visualizing results as summary"
#' subtitle: "MR of SH on GE"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Load files ####
#' ***
load("../results/02_PGS_04_PGSresults.RData")
PGS.result.genomewide = copy(PGS.result)
PGS.result.genomewide[,SNPsetting := "genome-wide"]

load("../results/02_PGS_06_PGSresults_pathway.RData")
PGS.result.pathwaywide = copy(PGS.result)
PGS.result.pathwaywide[,SNPsetting := "pathway-wide"]

load("../results/02_PGS_08_PGSresults_SERPINA6.RData")
PGS.result.SERPINA6 = copy(PGS.result)
PGS.result.SERPINA6[,SNPsetting := "SERPINA6"]

#' # Filter results ####
#' ***
#' - TT: men and women, p-value threshold 0.001, pathway & genome-wide
#' - E2: men and women, p-value threshold 0.001, pathway & genome-wide
#' - CORT: combined, p-value threshold 0.001 & 1e-4, SERPINA
#' 
PGS.result1 = rbind(PGS.result.genomewide,PGS.result.pathwaywide)
PGS.result1 = PGS.result1[range=="001",]

PGS.result2 = copy(PGS.result.SERPINA6)
PGS.result2 = PGS.result2[range %in% c("1e-3","1e-4")]

PGS.result = rbind(PGS.result1,PGS.result2)

#' # Create plots ####
#' ***
#' I want five plots, one per hormone and setting. Adult and Heart, and SNP setting is to be plotted in one plot. 
#' 
PGS.result[,trait := paste(hormone,setting,sep="_")]
myTraits = unique(PGS.result$trait)

PGS.result[,printP := signif(pval, digits = 3)]
PGS.result[,printP2 := sub("e", "*x*10^", printP)]

PGS.result[,printRange := paste(study,SNPsetting, sep=" - ")]
PGS.result[,printRange := gsub("genome-wide","gw",printRange)]
PGS.result[,printRange := gsub("pathway-wide","pw",printRange)]

PGS.result[hormone == "CORT",printRange := paste(study,range, sep=" - ")]

for(i in 1:length(myTraits)){
  #i=1
  plotData = copy(PGS.result)
  plotData = PGS.result[trait == myTraits[i]]
  
  myMax = plotData[,max(R2)]
  myXaxis = "Study - SNP selection"
  if(myTraits[i]=="CORT_combined") myXaxis = "Study - p-value threshold"
  
  ggplot(data = plotData, aes(x = factor(printRange), y = R2)) +
    # Specify that we want to print p-value on top of the bars
    geom_text(
      aes(label = paste(printP2)),
      vjust = -1.5,
      hjust = 0,
      angle = 45,
      cex = 4,
      parse = T
    )  +
    # Specify the range of the plot, *1.25 to provide enough space for the p-values
    scale_y_continuous(limits = c(0, myMax * 1.25)) +
    # Specify the axis labels
    #xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    xlab(myXaxis) +
    ylab(expression(paste("PGS model fit:  ", R ^ 2))) +
    # Draw a bar plot
    geom_bar(aes(fill = -log10(pval)), stat = "identity") +
    # Specify the colors
    scale_fill_gradient2(
      low = "dodgerblue",
      high = "firebrick",
      mid = "dodgerblue",
      midpoint = 1e-4,
      name = bquote(atop(-log[10] ~ model, italic(P) - value),)
    ) +
    # Some beautification of the plot
    theme_classic() + theme(
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  # save the plot
  ggsave(paste0("../results/02_PGS_09_PGS_barplots/",myTraits[i],".png"), height = 7, width = 7)
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
