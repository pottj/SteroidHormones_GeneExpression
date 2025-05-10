#' ---
#' title: "TWAS - Visualizing results (bar plot and vulcano plots)"
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
#' In LIFE-Adult, I want to visualize the  
#' 
#' - number of associated genes per models (BMI adjustment vs. no adjustment)
#' - foldchange by p-value (vulcano plot, per trait) 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/03_TWAS_01_SummaryStatistics_LIFEAdult_hierFDR.RData")
Adult = copy(myTab)

#' # Bar plots ####
#' *** 
#' 
#' - Plot 1: no BMI adjustment
#' - Plot 2: BMI adjustment
#' - Color coding per hormone and setting 
#'    - CORT: green 
#'    - TT: blue 
#'    - E2: red
#'    - women: light
#'    - men: dark
#'    
#' Number of associated genes per trait and setting
plotData = Adult[hierFDR==T,.N,by=c("phenotype","setting","model")]
plotData
plotData[,ratio := N/20972]
plotData[,phenotype := gsub("TESTO","TT",phenotype)]
plotData[,setting := gsub("combined","comb",setting)]
plotData[,trait2 := paste(phenotype,setting,model,sep=" - ")]
plotData[,trait := paste(phenotype,setting,sep=" - ")]

myMax = plotData[,max(N)]

ggplot(data = plotData[model=="noBMIadj"], aes(x = trait, y = N, fill = trait)) +
  # Draw a bar plot
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  # Specify that we want to print N on top of the bars
  geom_text(
    aes(label = paste(N)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, myMax * 1.25)) +
  # Specify the colors
  scale_fill_manual(values = c("#47D45A","#F2AA84","#FBE3D6","#61CBF4","#CAEEFB"), name = "Hormones") +
  # Specify the axis labels
  xlab(expression("Hormone - sample set")) +
  ylab(expression("Associated genes (FDR)")) +
  # Some beautification of the plot
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
# save the plot
ggsave(paste0("../results/03_TWAS_04_Plots/Barplot_noBMIadj.png"), height = 7, width = 7)

ggplot(data = plotData[model=="BMIadj"], aes(x = trait, y = N, fill = trait)) +
  # Draw a bar plot
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  # Specify that we want to print N on top of the bars
  geom_text(
    aes(label = paste(N)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, myMax * 1.25)) +
  # Specify the colors
  scale_fill_manual(values = c("#47D45A","#F2AA84","#FBE3D6","#61CBF4","#CAEEFB"), name = "Hormones") +
  # Specify the axis labels
  xlab(expression("Hormone - sample set")) +
  ylab(expression("Associated genes (FDR)")) +
  # Some beautification of the plot
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
# save the plot
ggsave(paste0("../results/03_TWAS_04_Plots/Barplot_BMIadj.png"), height = 7, width = 7)

#' # Volcano plots ####
#' ***
#' 
Adult[,diffexpressed := "NO"]
Adult[beta > 0.05 & hierFDR==T,diffexpressed := "UP"]
Adult[beta < -0.05 & hierFDR==T,diffexpressed := "DOWN"]

Adult[,label := ""]
Adult[diffexpressed != "NO",label := symbol_INGENUITY]

myTraits = unique(Adult$trait)
Adult[,trait2 := paste(phenotype,setting,sep=" - ")]

ggplot(data=Adult[phenotype == "CORT",],
       aes(x=beta, y=-log10(P.Value.adj2), col=diffexpressed, label=label)) +
  facet_grid(trait2 ~ model, axes = "all", axis.labels = "all_x",scales = "free")+
  geom_point() +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.05, 0.05), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")+
  xlab(expression("log2(Foldchange)")) +
  ylab(expression("-log10(adjusted p-value)"))+
  theme(legend.position = "none")

ggsave(paste0("../results/03_TWAS_04_Plots/VolcanoPlot_CORT.png"), height = 7, width = 14)

ggplot(data=Adult[phenotype == "TESTO",],
       aes(x=beta, y=-log10(P.Value.adj2), col=diffexpressed, label=label)) +
  facet_grid(trait2 ~ model, axes = "all", axis.labels = "all_x",scales = "free")+
  geom_point() +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.05, 0.05), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")+
  xlab(expression("log2(Foldchange)")) +
  ylab(expression("-log10(adjusted p-value)"))+
  theme(legend.position = "none")

ggsave(paste0("../results/03_TWAS_04_Plots/VolcanoPlot_TESTO.png"), height = 14, width = 14)

ggplot(data=Adult[phenotype == "E2",],
       aes(x=beta, y=-log10(P.Value.adj2), col=diffexpressed, label=label)) +
  facet_grid(trait2 ~ model, axes = "all", axis.labels = "all_x",scales = "free")+
  geom_point() +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.05, 0.05), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")+
  xlab(expression("log2(Foldchange)")) +
  ylab(expression("-log10(adjusted p-value)"))+
  theme(legend.position = "none")

ggsave(paste0("../results/03_TWAS_04_Plots/VolcanoPlot_E2.png"), height = 14, width = 14)

for(i in 1:length(myTraits)){
  #i=1
  ggplot(data=Adult[trait == myTraits[i],], 
         aes(x=beta, y=-log10(P.Value.adj2), col=diffexpressed, label=label)) +
    geom_point() + 
    theme_classic() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.05, 0.05), col="grey") +
    geom_hline(yintercept=-log10(0.05), col="grey")+
    xlab(expression("log2(Foldchange)")) + 
    ylab(expression("-log10(adjusted p-value)"))+ 
    theme(legend.position = "none")
  ggsave(paste0("../results/03_TWAS_04_Plots/VolcanoPlot_",myTraits[i],".png"), height = 7, width = 7)
  
}

#' 
#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
