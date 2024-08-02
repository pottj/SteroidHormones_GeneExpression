#' ---
#' title: "Get scores"
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
#' In this script, I will create the PLINK2 calls to generate the PGS of TT, E2 and CORT using the weights from the publicly available GWAS data. 
#' 
#' As sanity check, I test the association between PGS and hormone data. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_angmar.R")
.libPaths()

#' # Scores Adult ####
#' ***
wheights = list.files(path = "../temp/", pattern = "ScoreWeights")
wheights

load(paste0(path_LIFEprepped,"01_LIFEAdult_filtered_final.RData"))
myTab

dumTab = foreach(i = 1:length(wheights))%do%{
  #i=1
  dummy1 = gsub("_ScoreWeights.*","",wheights[i])
  dummy1 = gsub(".*_","",dummy1)
  
  dummy2 = gsub(".*ScoreWeights_","",wheights[i])
  dummy2 = gsub(".txt","",dummy2)
  
  myCall1 = paste0(path_plink2,
                   " --pfile ",path_LIFEAdult_Genetics,
                   " --score ../temp/",wheights[i],
                   " --out ../results/05_Scores_Adult_",dummy1,"_",dummy2)
  myCall1
  system(myCall1)
  
  myScore = fread(paste0("../results/05_Scores_Adult_",dummy1,"_",dummy2,".sscore"))
  matched = match(myTab$ALIQUOT_genetics,myScore$IID)
  table(is.na(matched))
  myTab[,score := myScore[matched,SCORE1_AVG]]
  setnames(myTab,"score",paste0("PGS_",dummy1,"_",dummy2))
  
}

# some checks
myTab[,PGS_TT_best:=PGS_TT_Men_best]
myTab[GENDER==2,PGS_TT_best:=PGS_TT_Women_best]
myTab[,PGS_TT_same:=PGS_TT_Men_same]
myTab[GENDER==2,PGS_TT_same:=PGS_TT_Women_same]
summary(lm(TESTO ~ PGS_TT_Men_best, data=myTab,subset=GENDER==1))
summary(lm(TESTO ~ PGS_TT_Men_same, data=myTab,subset=GENDER==1))
summary(lm(TESTO ~ PGS_TT_Women_best, data=myTab,subset=GENDER==2))
summary(lm(TESTO ~ PGS_TT_Women_same, data=myTab,subset=GENDER==2))
summary(lm(log(TESTO) ~ PGS_TT_best + GENDER, data=myTab))
summary(lm(log(TESTO) ~ PGS_TT_same + GENDER, data=myTab,group!="premenopausal"))

myTab[,PGS_E2_best:=PGS_E2_Men_best]
myTab[GENDER==2,PGS_E2_best:=PGS_E2_Women_best]
myTab[,PGS_E2_same:=PGS_E2_Men_same]
myTab[GENDER==2,PGS_E2_same:=PGS_E2_Women_same]
summary(lm(log(E2) ~ PGS_E2_Men_best, data=myTab,subset=GENDER==1))
summary(lm(log(E2) ~ PGS_E2_Men_same, data=myTab,subset=GENDER==1))
summary(lm(log(E2) ~ PGS_E2_Women_best, data=myTab,subset=GENDER==2 & group=="postmenopausal"))
summary(lm(log(E2) ~ PGS_E2_Women_same, data=myTab,subset=GENDER==2 & group=="postmenopausal"))
summary(lm(log(E2) ~ PGS_E2_best + GENDER, data=myTab,group!="premenopausal"))
summary(lm(log(E2) ~ PGS_E2_same + GENDER, data=myTab,group!="premenopausal"))

summary(lm(log(CORT) ~ PGS_CORT_Combined_best + GENDER + D126_time + log(D074_BMI) + D141_smokeStatus, data=myTab))

save(myTab,file=paste0(path_LIFEprepped,"05_LIFEAdult_filtered_final_PGS.RData"))

#' # Scores Heart ####
#' ***
wheights = list.files(path = "../temp/", pattern = "ScoreWeights")
wheights

load(paste0(path_LIFEprepped,"02_LIFEHeart_filtered_final.RData"))
myTab

dumTab = foreach(i = 1:length(wheights))%do%{
  #i=1
  dummy1 = gsub("_ScoreWeights.*","",wheights[i])
  dummy1 = gsub(".*_","",dummy1)
  
  dummy2 = gsub(".*ScoreWeights_","",wheights[i])
  dummy2 = gsub(".txt","",dummy2)
  
  myCall1 = paste0(path_plink2,
                   " --pfile ",path_LIFEHeart_Genetics,
                   " --score ../temp/",wheights[i],
                   " --out ../results/05_Scores_Heart_",dummy1,"_",dummy2)
  myCall1
  system(myCall1)
  
  myScore = fread(paste0("../results/05_Scores_Heart_",dummy1,"_",dummy2,".sscore"))
  matched = match(myTab$ALIQUOT_genetics,myScore$`#IID`)
  table(is.na(matched))
  myTab[,score := myScore[matched, SCORE1_AVG]]
  setnames(myTab,"score",paste0("PGS_",dummy1,"_",dummy2))
  
}

# some checks
myTab[,PGS_TT_best:=PGS_TT_Men_best]
myTab[GENDER==2,PGS_TT_best:=PGS_TT_Women_best]
myTab[,PGS_TT_same:=PGS_TT_Men_same]
myTab[GENDER==2,PGS_TT_same:=PGS_TT_Women_same]
summary(lm(TESTO ~ PGS_TT_Men_best, data=myTab,subset=GENDER==1))
summary(lm(TESTO ~ PGS_TT_Men_same, data=myTab,subset=GENDER==1))
summary(lm(TESTO ~ PGS_TT_Women_best, data=myTab,subset=GENDER==2 & group!="premenopausal"))
summary(lm(TESTO ~ PGS_TT_Women_same, data=myTab,subset=GENDER==2 & group!="premenopausal"))
summary(lm(log(TESTO) ~ PGS_TT_best + GENDER, data=myTab,group!="premenopausal"))
summary(lm(log(TESTO) ~ PGS_TT_same + GENDER, data=myTab,group!="premenopausal"))

myTab[,PGS_E2_best:=PGS_E2_Men_best]
myTab[GENDER==2,PGS_E2_best:=PGS_E2_Women_best]
myTab[,PGS_E2_same:=PGS_E2_Men_same]
myTab[GENDER==2,PGS_E2_same:=PGS_E2_Women_same]
summary(lm(log(E2) ~ PGS_E2_Men_best, data=myTab,subset=GENDER==1))
summary(lm(log(E2) ~ PGS_E2_Men_same, data=myTab,subset=GENDER==1))
summary(lm(log(E2) ~ PGS_E2_Women_best, data=myTab,subset=GENDER==2 & group=="postmenopausal"))
summary(lm(log(E2) ~ PGS_E2_Women_same, data=myTab,subset=GENDER==2 & group=="postmenopausal"))
summary(lm(log(E2) ~ PGS_E2_best + GENDER, data=myTab,group!="premenopausal"))
summary(lm(log(E2) ~ PGS_E2_same + GENDER, data=myTab,group!="premenopausal"))

summary(lm(log(CORT) ~ PGS_CORT_Combined_best + T991_time + log(D157_BMI), data=myTab))

save(myTab,file=paste0(path_LIFEprepped,"05_LIFEHeart_filtered_final_PGS.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
