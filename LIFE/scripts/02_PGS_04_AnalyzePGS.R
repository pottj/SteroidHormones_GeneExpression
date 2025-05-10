#' ---
#' title: "PGS - Analysing PGS and Visualizing results"
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
#' I want to load all the scores per study and create one score file. As template I use the respective "samplesAsImputed" file, as it also contains the genetic PCs. 
#' 
#' Then I load the phenotype files, and perform linear regression of the scores on the respective outcome in the respective subset of samples, using either no adjustment, or adjusting for age, BMI, smoking, and the PCs. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Create input parameters ####
#' ***
files = list.files(path = path_SumStats_QC)
files = files[grepl(".gz",files)]
files = files[!grepl("PathwayRegions",files)]
files = files[!grepl("SERPINA6",files)]
studies = c("Adult","Heart")

#' # Loop to load scores ####
#' ***
for(j in 1:length(studies)){
  #j=1
  message("Working on study LIFE-",studies[j])
  
  # Step 1: load samplesAsImputed
  if(studies[j]=="Adult"){
    loaded = load(path_LIFEAdult_GeneticPCs)
    SamplesAsImputed = get(loaded)
    SamplesAsImputed = SamplesAsImputed[filterGonosom==F,]
    SamplesAsImputed = SamplesAsImputed[,c(2,5,8:17)]
    #names(SamplesAsImputed)
  }else{
    loaded = load(path_LIFEHeart_GeneticPCs)
    SamplesAsImputed = get(loaded)
    setDT(SamplesAsImputed)
    SamplesAsImputed = SamplesAsImputed[filterGonosom==F,]
    SamplesAsImputed = SamplesAsImputed[,c(6,10,12:21)]
    setnames(SamplesAsImputed,"Aliquot_Genetik","Aliquot")
    #names(SamplesAsImputed)
  }
  
  # Step 2: get score files per study
  myScoreFiles = list.files(path="../results/02_PGS_03_PLINK_PGS/", pattern="sscore")
  myScoreFiles = myScoreFiles[grepl(studies[j],myScoreFiles)]
  
  # Step 3: loop over the 35 scores
  for(i in 1:length(myScoreFiles)){
    #i=1
    dummy = myScoreFiles[i]
    dummy = gsub(studies[j],"",dummy)
    dummy = gsub("_Score_","",dummy)
    dummy = gsub(".sscore","",dummy)
    hormone = gsub("_.*","",dummy)
    setting = gsub("[.].*","",dummy)
    setting = gsub(".*_","",setting)
    pval_range = gsub(".*[.]","",dummy)
    message("       and hormone ",hormone," in setting ",setting," with max. p-value of 0.",pval_range)
    
    # Step 4: load score
    myScore = fread(paste0("../results/02_PGS_03_PLINK_PGS/",myScoreFiles[i]))
    
    # Step 5: add score to samplesAsImputed
    matched = match(SamplesAsImputed$Aliquot,myScore$IID) 
    table(SamplesAsImputed$Aliquot == myScore$IID[matched])
    SamplesAsImputed[,score := myScore[matched,SCORE1_AVG]]
    
    # Step 6: rename score column: hormone_setting_range
    setnames(SamplesAsImputed,"score",paste(hormone,setting,pval_range,sep="_"))
    
  }
  
  # Step 7: save object 
  save(SamplesAsImputed, file = paste0(path_LIFEprepped,"genetics/",studies[j],"_PGS.RData"))
  
}

load(paste0(path_LIFEprepped,"genetics/Adult_PGS.RData"))
myTab_Adult_PGS = copy(SamplesAsImputed)

load(paste0(path_LIFEprepped,"genetics/Heart_PGS.RData"))
myTab_Heart_PGS = copy(SamplesAsImputed)

#' # Merge with phenotype data ####
#' ***
load(paste0(path_LIFEprepped,"phenotypes/Adult_QC.RData"))
myTab_Adult = copy(myTab) 

load(paste0(path_LIFEprepped,"phenotypes/Heart_QC.RData"))
myTab_Heart = copy(myTab) 

#' I merge the scores to the phenotype file. **I remove the score for CORT_combined_001 as it is always NA.**
matched = match(myTab_Adult$ALIQUOT_genetics,myTab_Adult_PGS$Aliquot)
table(myTab_Adult$ALIQUOT_genetics == myTab_Adult_PGS$Aliquot[matched])
table(myTab_Adult$GENDER,myTab_Adult_PGS$sex[matched])
myTab_Adult = cbind(myTab_Adult,myTab_Adult_PGS[matched,-c(1,2,13),with=F])

matched = match(myTab_Heart$ALIQUOT_genetics,myTab_Heart_PGS$Aliquot)
table(myTab_Heart$ALIQUOT_genetics == myTab_Heart_PGS$Aliquot[matched])
table(myTab_Heart$GENDER,myTab_Heart_PGS$sex[matched])
myTab_Heart = cbind(myTab_Heart,myTab_Heart_PGS[matched,-c(1,2,13),with=F])

#' # Run linear regression ####
#' ***
myScores = names(myTab_Adult)[42:75]

dumTab1 = foreach(i = 1:length(myScores))%do%{
  #i=1
  
  # Step 1: get hormone and setting out of score name
  hormone = gsub("_.*","",myScores[i])
  range = gsub(".*_","",myScores[i])
  setting = gsub(hormone,"",myScores[i])
  setting = gsub(range,"",setting)
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
                      range = range,
                      R2 = prs.r2,
                      sampleSize = length(model2$residuals),
                      beta = prs.coef[1],
                      SE = prs.coef[2],
                      tval = prs.coef[3],
                      pval = prs.coef[4])
  result
}
PGS.result.Adult = rbindlist(dumTab1)
PGS.result.Adult[,max(R2),by=c("hormone","setting")]
PGS.result.Adult[,min(pval),by=c("hormone","setting")]
PGS.result.Adult[pval<5e-8,]

#' Only strong instruments for testosterone. 
#' 
#' Now the same for Heart (some column names are different)
#' 
dumTab2 = foreach(i = 1:length(myScores))%do%{
  #i=1
  
  # Step 1: get hormone and setting out of score name
  hormone = gsub("_.*","",myScores[i])
  range = gsub(".*_","",myScores[i])
  setting = gsub(hormone,"",myScores[i])
  setting = gsub(range,"",setting)
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
                      range = range,
                      R2 = prs.r2,
                      sampleSize = length(model2$residuals),
                      beta = prs.coef[1],
                      SE = prs.coef[2],
                      tval = prs.coef[3],
                      pval = prs.coef[4])
  result
}
PGS.result.Heart = rbindlist(dumTab2)
PGS.result.Heart[,max(R2),by=c("hormone","setting")]
PGS.result.Heart[,min(pval),by=c("hormone","setting")]
PGS.result.Heart[pval<5e-8,]
PGS.result.Heart[hormone=="TT"]

#' In the Heart study only testosterone in women has strong instruments. 
#' 
#' # Save PGS results ####
#' ***
PGS.result.Adult[,study := "Adult"]
PGS.result.Heart[,study := "Heart"]

PGS.result = rbind(PGS.result.Adult, PGS.result.Heart)

save(PGS.result,file = "../results/02_PGS_04_PGSresults.RData")
write.table(PGS.result, file = "../results/02_PGS_04_PGSresults.txt",quote = F,sep="\t",row.names = F,col.names = T,dec = ".")

#' # Make some plots ####
#' ***
#' For all hormones and settings. First, I generate a pretty format for p-value output
#'  
plotdata = copy(PGS.result)

plotdata[,printP := signif(pval, digits = 3)]
plotdata[,printP2 := sub("e", "*x*10^", printP)]
plotdata[,printRange := paste0("0.",range)]

#' Now I need a loop over each trait 
ToDoList = plotdata[,.N,by=c("hormone","setting","study")]

for(i in 1:dim(ToDoList)[1]){
  #i=1
  myMax = plotdata[hormone == ToDoList[i,hormone] & setting == ToDoList[i,setting],max(R2)]
  ggplot(data = plotdata[hormone == ToDoList[i,hormone] & setting == ToDoList[i,setting] & study == ToDoList[i,study]], 
         aes(x = factor(printRange), y = R2)) +
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
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
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
  filename = paste(ToDoList[i,hormone],ToDoList[i,setting],ToDoList[i,study],sep="_")
  ggsave(paste0("../results/02_PGS_04_PGS_barplots/",filename,".png"), height = 7, width = 7)
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
