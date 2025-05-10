#' ---
#' title: "PGS - Selecting the scores"
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
load(paste0(path_LIFEprepped,"genetics/Adult_PGS.RData"))
myTab_Adult_PGS_genomewide = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Heart_PGS.RData"))
myTab_Heart_PGS_genomewide = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Adult_PGS_pathway.RData"))
myTab_Adult_PGS_pathway = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Heart_PGS_pathway.RData"))
myTab_Heart_PGS_pathway = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Adult_PGS_SERPINA6.RData"))
myTab_Adult_PGS_SERPINA6 = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Heart_PGS_SERPINA6.RData"))
myTab_Heart_PGS_SERPINA6 = copy(SamplesAsImputed)

#' # Merge and filter data ####
#' ***
myTab_Adult_PGS = cbind(myTab_Adult_PGS_genomewide[,c(1:12,34,41),with=F],myTab_Adult_PGS_pathway[,20],myTab_Adult_PGS_SERPINA6[,13])
names(myTab_Adult_PGS)[c(13:16)] = c("TT_men","TT_women","E2_men","CORT_comb")

myTab_Heart_PGS = cbind(myTab_Heart_PGS_genomewide[,c(1:12,34,41),with=F],myTab_Heart_PGS_pathway[,20],myTab_Heart_PGS_SERPINA6[,13])
names(myTab_Heart_PGS)[c(13:16)] = c("TT_men","TT_women","E2_men","CORT_comb")

save(myTab_Adult_PGS, file = paste0(path_LIFEprepped,"genetics/Adult_PGS_selection.RData"))
save(myTab_Heart_PGS, file = paste0(path_LIFEprepped,"genetics/Heart_PGS_selection.RData"))

#' # Rerun association with scaled scores ####
#' ***
myScores = names(myTab_Adult_PGS)[13:16]

#' ## Adult ####
load(paste0(path_LIFEprepped,"phenotypes/Adult_QC.RData"))
myTab_Adult = copy(myTab) 
matched = match(myTab_Adult$ALIQUOT_genetics,myTab_Adult_PGS$Aliquot)
table(myTab_Adult$ALIQUOT_genetics == myTab_Adult_PGS$Aliquot[matched])
table(myTab_Adult$GENDER,myTab_Adult_PGS$sex[matched])
myTab_Adult = cbind(myTab_Adult,myTab_Adult_PGS[matched,-c(1,2),with=F])

dumTab1 = foreach(i = 1:length(myScores))%do%{
  #i=1
  
  # Step 1: get hormone and setting out of score name
  hormone = gsub("_.*","",myScores[i])
  setting = gsub(hormone,"",myScores[i])
  setting = gsub("_","",setting)
  
  # Step 2: filter for correct setting
  data = copy(myTab_Adult)
  if(setting == "men"){
    data = data[group=="men",]
  }else if(setting == "women"){
    data = data[group!="men",]
  }
  
  # Step 3: identify hormone and score
  if(hormone == "CORT"){
    data[,hormone := CORT]
  }else if(hormone == "TT"){
    data[,hormone := TESTO]
  }else if(hormone == "E2"){
    data[,hormone := E2]
  }
  data[,score := get(myScores[i])]
  data[,score := scale(score)]
  
  # Step 4: do linear regression
  model0 = lm(log(hormone) ~ GENDER + AGE + D126_time + log(D074_BMI), data = data)
  if(setting != "men"){
    model2 = lm(log(hormone) ~ score + GENDER + AGE + group + D126_time + log(D074_BMI) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data)
  }else{
    model2 = lm(log(hormone) ~ score + GENDER + AGE + D126_time + log(D074_BMI) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data)
  }
  
  # Step 5: summarize results
  null.r2 = summary(model0)$r.squared
  model.r2 = summary(model2)$r.squared
  prs.r2 <- model.r2-null.r2
  prs.coef <- summary(model2)$coeff["score",]
  
  result = data.table(score = myScores[i],
                      hormone = hormone, 
                      setting = setting,
                      range = 0.001,
                      R2 = prs.r2,
                      sampleSize = length(model2$residuals),
                      beta = prs.coef[1],
                      SE = prs.coef[2],
                      tval = prs.coef[3],
                      pval = prs.coef[4])
  result
}
PGS.result.Adult = rbindlist(dumTab1)
PGS.result.Adult

#' ## Heart ####
load(paste0(path_LIFEprepped,"phenotypes/Heart_QC.RData"))
myTab_Heart = copy(myTab) 
matched = match(myTab_Heart$ALIQUOT_genetics,myTab_Heart_PGS$Aliquot)
table(myTab_Heart$ALIQUOT_genetics == myTab_Heart_PGS$Aliquot[matched])
table(myTab_Heart$GENDER,myTab_Heart_PGS$sex[matched])
myTab_Heart = cbind(myTab_Heart,myTab_Heart_PGS[matched,-c(1,2),with=F])

dumTab2 = foreach(i = 1:length(myScores))%do%{
  #i=1
  
  # Step 1: get hormone and setting out of score name
  hormone = gsub("_.*","",myScores[i])
  range = gsub(".*_","",myScores[i])
  setting = gsub(hormone,"",myScores[i])
  setting = gsub("_","",setting)
  
  # Step 2: filter for correct setting
  data = copy(myTab_Heart)
  if(setting == "men"){
    data = data[group=="men",]
  }else if(setting == "women"){
    data = data[group!="men",]
  }
  
  # Step 3: identify hormone and score
  if(hormone == "CORT"){
    data[,hormone := CORT]
  }else if(hormone == "TT"){
    data[,hormone := TESTO]
  }else if(hormone == "E2"){
    data[,hormone := E2]
  }
  data[,score := get(myScores[i])]
  data[,score := scale(score)]
  
  # Step 4: do linear regression
  model0 = lm(log(hormone) ~ GENDER + AGE + T991_time + log(D157_BMI), data = data)
  if(setting != "men"){
    model2 = lm(log(hormone) ~ score + GENDER + AGE + group + T991_time + log(D157_BMI) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data)
  }else{
    model2 = lm(log(hormone) ~ score + GENDER + AGE + T991_time + log(D157_BMI) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data)
  }
  
  # Step 5: summarize results
  null.r2 = summary(model0)$r.squared
  model.r2 = summary(model2)$r.squared
  prs.r2 <- model.r2-null.r2
  prs.coef <- summary(model2)$coeff["score",]
  
  result = data.table(score = myScores[i],
                      hormone = hormone, 
                      setting = setting,
                      range = 0.001,
                      R2 = prs.r2,
                      sampleSize = length(model2$residuals),
                      beta = prs.coef[1],
                      SE = prs.coef[2],
                      tval = prs.coef[3],
                      pval = prs.coef[4])
  result
}
PGS.result.Heart = rbindlist(dumTab2)
PGS.result.Heart

#' # Save PGS results ####
#' ***
PGS.result.Adult[,study := "Adult"]
PGS.result.Heart[,study := "Heart"]

PGS.result = rbind(PGS.result.Adult, PGS.result.Heart)

save(PGS.result,file = "../results/02_PGS_10_PGSresults_selectionScaled.RData")
fwrite(PGS.result, file = "../results/02_PGS_10_PGSresults_selectionScaled.txt",quote = F,sep="\t",row.names = F,col.names = T,dec = ".")

#' # Make some plots ####
#' ***
#' For all hormones and settings. First, I generate a pretty format for p-value output
#'  
plotdata = copy(PGS.result)

plotdata[,printP := signif(pval, digits = 3)]
plotdata[,printP2 := sub("e", "*x*10^", printP)]
plotdata[,trait := gsub("_", " - ", score)]

myMax = plotdata[study == "Adult",max(R2)]
ggplot(data = plotdata[study == "Adult"], 
       aes(x = factor(trait), y = R2)) +
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
  xlab("Horomone - sample setting") +
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
ggsave(paste0("../results/02_PGS_09_PGS_barplots/Adult_finalScores.png"), height = 7, width = 7)

myMax = plotdata[study == "Heart",max(R2)]
ggplot(data = plotdata[study == "Heart"], 
       aes(x = factor(trait), y = R2)) +
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
  xlab("Horomone - sample setting") +
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
ggsave(paste0("../results/02_PGS_09_PGS_barplots/Heart_finalScores.png"), height = 7, width = 7)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

