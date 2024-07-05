#' ---
#' title: "Run limma in LIFE-Heart"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_forostar.R")
.libPaths()

source(func_TWAS)
source(func_TWASExtract)

#' # Load data ####
#' ***
#' Load the phenotype data (hormones and covariates) and the Expression Set (GE)
#' 
load(paste0(path_LIFEprepped,"02_LIFEHeart_filtered_final.RData"))
myTab

loaded1 = load(path_LIFEHeart_GE)
loaded1
eset_B3

loaded2 = load(path_LIFEAdult_GE)
loaded2
eset_A1

#' # Transform data ####
#' ***
#' I want to log-transform the hormone data and BMI, make smoking binary (active smoking only), and change some column names
#' 
names(myTab)
myTab[,CORT := log(CORT)]
myTab[,TESTO := log(TESTO)]
myTab[,E2 := log(E2)]
myTab[,D157_BMI := log(D157_BMI)]

myTab[,table(D222_smokeStatus)]
myTab[D222_smokeStatus==1,D222_smokeStatus:=0]
myTab[D222_smokeStatus==2,D222_smokeStatus:=1]

setnames(myTab,"GENDER","SEX")
setnames(myTab,"T991_time","TimeBloodSample")
setnames(myTab,"T991_fasting","FastingHours")
setnames(myTab,"D157_BMI","BMI")
setnames(myTab,"D222_smokeStatus","activeSmoking")

myTab = myTab[TWAS==T,]

#' # Merge data ####
#' ***
#' Add the relevant columns to the Expression Set
#' 
temp1 = pData(eset_B3)

toadd = c("CORT","TESTO","E2" ,"activeSmoking","BMI","group","TimeBloodSample","OUremoved")
table(is.element(myTab$ALIQUOT_GE,temp1$Aliquot))
table(is.element(temp1$Aliquot,myTab$ALIQUOT_GE))

for(i in toadd) {
  pData(eset_B3)[i]    =  myTab[ match_hk(pData(eset_B3)$Aliquot,  myTab$ALIQUOT_GE), i, with = F]
}

#' # Filter data ####
#' ***
#' Reduce Expression Set to samples with hormone data 
eset_SH = eset_B3[,!is.na(eset_B3$CORT)]
dim(eset_SH)
temp2 = pData(eset_SH)
table(is.na(temp2$CORT),temp2$group)

#' Select high quality probes only (same as in Adult)
temp3 = fData(eset_A1)
filt_features<-temp3$probe_QCok==T & temp3$unique_mapper==T
table(filt_features)
goodfeatures<-temp3$PROBE_ID[filt_features]

temp4 = fData(eset_SH)
filt_features<- is.element(temp4$PROBE_ID,goodfeatures)
table(filt_features)
goodfeatures<-temp4$PROBE_ID[filt_features]
eset_SH = eset_SH[goodfeatures,]
dim(eset_SH)

#' save filtered eSet
#' 
save(eset_SH, file = paste0(path_LIFEprepped,"eSet_Heart_TWAS.RData"))

#' # Run limma #### 
#' ***
SHs = c("CORT","TESTO","E2")

#' ## Men ####
gx_assoc_men_filt = calculateLimma(todovars = SHs, 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

gx_assoc_men_filt2= calculateLimma(todovars = SHs, 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

#' ## Postmenopausal women ####
gx_assoc_postWomen_filt = calculateLimma(todovars = SHs, 
                                         conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                         eset = eset_SH[,eset_SH$group=="postmenopausal"], 
                                         doscale = T) 

gx_assoc_postWomen_filt2= calculateLimma(todovars = SHs, 
                                         conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                         eset = eset_SH[,eset_SH$group=="postmenopausal"], 
                                         doscale = T) 

#' ## Premenopausal women ####
gx_assoc_preWomen_filt = calculateLimma(todovars = SHs, 
                                        conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                        eset = eset_SH[,eset_SH$group=="premenopausal"], 
                                        doscale = T) 

gx_assoc_preWomen_filt2= calculateLimma(todovars = SHs, 
                                        conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                        eset = eset_SH[,eset_SH$group=="premenopausal"], 
                                        doscale = T) 

#' ## Women ####
gx_assoc_Women_filt = calculateLimma(todovars = SHs, 
                                     conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group"), 
                                     eset = eset_SH[,eset_SH$group!="men"], 
                                     doscale = T) 

gx_assoc_Women_filt2= calculateLimma(todovars = SHs, 
                                     conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group","BMI"), 
                                     eset = eset_SH[,eset_SH$group!="men"], 
                                     doscale = T) 

#' ## All ####
gx_assoc_all_filt = calculateLimma(todovars = SHs, 
                                     conf = c("SEX","AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                     eset = eset_SH, 
                                     doscale = T) 

gx_assoc_all_filt2= calculateLimma(todovars = SHs, 
                                     conf = c("SEX","AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                     eset = eset_SH, 
                                     doscale = T) 

#' # Save lists ####
#' ***
save(gx_assoc_men_filt, file="../results/02_TWAS_Heart_men.RData")
save(gx_assoc_men_filt2, file="../results/02_TWAS_Heart_men_BMIadj.RData")

save(gx_assoc_postWomen_filt, file="../results/02_TWAS_Heart_postWomen.RData")
save(gx_assoc_postWomen_filt2, file="../results/02_TWAS_Heart_postWomen_BMIadj.RData")

save(gx_assoc_preWomen_filt, file="../results/02_TWAS_Heart_preWomen.RData")
save(gx_assoc_preWomen_filt2, file="../results/02_TWAS_Heart_preWomen_BMIadj.RData")

save(gx_assoc_Women_filt, file="../results/02_TWAS_Heart_women.RData")
save(gx_assoc_Women_filt2, file="../results/02_TWAS_Heart_women_BMIadj.RData")

save(gx_assoc_all_filt, file="../results/02_TWAS_Heart_all.RData")
save(gx_assoc_all_filt2, file="../results/02_TWAS_Heart_all_BMIadj.RData")

#' # Extract Summary statistics
#' ***
#' I do not like lists - I want a data table!
#' 
#' ## No BMI adjustment ####

dumTab = foreach(i=1:length(SHs))%do%{
  #i=1
  message("Working on i=",i,", steroid hormone ",SHs[i])
  mySH = SHs[i]
  tab_1 = myExtractionFunction(gx_assoc_preWomen_filt,mySH,"women_pre")
  tab_2 = myExtractionFunction(gx_assoc_postWomen_filt,mySH,"women_post")
  tab_3 = myExtractionFunction(gx_assoc_Women_filt,mySH,"women")
  tab_4 = myExtractionFunction(gx_assoc_men_filt,mySH,"men")
  tab_5 = myExtractionFunction(gx_assoc_all_filt,mySH,"all")
  tab_heart = rbind(tab_1,tab_2,tab_3,tab_4,tab_5)
  tab_heart
}

Heart = rbindlist(dumTab)  
Heart[,study:="Heart"]
Heart[,uniqueID := paste(phenotype,setting,PROBE_ID,sep="::")]
table(duplicated(Heart$uniqueID))

names(Heart)

Heart = Heart[,c(18,22,23,24,1,2,14,20,13,19,21,15:17,9:11,3:8,12)]
Heart

save(Heart, file = "../results/02_TWAS_SummaryStatistics_LIFEHeart.RData")

outfn = paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart.txt")
fwrite(Heart, file = outfn,quote = F,sep="\t")
R.utils::gzip(paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart.txt"), 
              destname = paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart.txt.gz"))

#' ## BMI adjustment ####

dumTab = foreach(i=1:length(SHs))%do%{
  #i=1
  message("Working on i=",i,", steroid hormone ",SHs[i])
  mySH = SHs[i]
  tab_1 = myExtractionFunction(gx_assoc_preWomen_filt,mySH,"women_pre")
  tab_2 = myExtractionFunction(gx_assoc_postWomen_filt,mySH,"women_post")
  tab_3 = myExtractionFunction(gx_assoc_Women_filt2,mySH,"women")
  tab_4 = myExtractionFunction(gx_assoc_men_filt2,mySH,"men")
  tab_5 = myExtractionFunction(gx_assoc_all_filt2,mySH,"all")
  tab_heart = rbind(tab_1,tab_2,tab_3,tab_4,tab_5)
  tab_heart
}

Heart2 = rbindlist(dumTab)  
Heart2[,study:="Heart_BMIadj"]
Heart2[,uniqueID := paste(phenotype,setting,PROBE_ID,sep="::")]
table(duplicated(Heart2$uniqueID))

names(Heart2)

Heart2 = Heart2[,c(18,22,23,24,1,2,14,20,13,19,21,15:17,9:11,3:8,12)]
Heart2

save(Heart2, file = "../results/02_TWAS_SummaryStatistics_LIFEHeart_BMIadj.RData")

outfn = paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart_BMIadj.txt")
fwrite(Heart2, file = outfn,quote = F,sep="\t")
R.utils::gzip(paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart_BMIadj.txt"), 
              destname = paste0("../results/02_TWAS_SummaryStatistics_LIFEHeart_BMIadj.txt.gz"))

#' # Hierarchical FDR ####
#' ***
myTab = copy(Heart)
myTab[,table(setting,phenotype)]

mySettings = unique(myTab$setting)

dumTab = foreach(i = 1:length(mySettings))%do%{
  #i=1
  mySetting = mySettings[i]
  myTab2 = copy(myTab)
  myTab2 = myTab2[setting == mySetting,]
  
  myFDR = addHierarchFDR(pvalues = myTab2$P.Value, categs = myTab2$phenotype)
  k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
  n = length(unique(myFDR[,category]))
  table(myFDR$hierarch_fdr5proz)
  myTab2[,P.Value.adj1 := myFDR$fdr_level1]
  myTab2[,hierFDR := myFDR$hierarch_fdr5proz]
  myTab2[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
  myTab2
}

myTab2 = rbindlist(dumTab)
myTab2[,table(hierFDR,P.Value.adj1<0.05)]
myTab2[,table(hierFDR,P.Value.adj2<0.05)]
myTab2[,table(hierFDR,setting)]
myTab2[hierFDR==T,table(phenotype,setting)]

save(myTab2, file = "../results/02_TWAS_SummaryStatistics_LIFEHeart_hierFDR.RData")

#' ## BMI adjusted ####
myTab = copy(Heart2)
myTab[,table(setting,phenotype)]

mySettings = unique(myTab$setting)

dumTab = foreach(i = 1:length(mySettings))%do%{
  #i=1
  mySetting = mySettings[i]
  myTab2 = copy(myTab)
  myTab2 = myTab2[setting == mySetting,]
  
  myFDR = addHierarchFDR(pvalues = myTab2$P.Value, categs = myTab2$phenotype)
  k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
  n = length(unique(myFDR[,category]))
  table(myFDR$hierarch_fdr5proz)
  myTab2[,P.Value.adj1 := myFDR$fdr_level1]
  myTab2[,hierFDR := myFDR$hierarch_fdr5proz]
  myTab2[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
  myTab2
}

myTab3 = rbindlist(dumTab)
myTab3[,table(hierFDR,P.Value.adj1<0.05)]
myTab3[,table(hierFDR,P.Value.adj2<0.05)]
myTab3[,table(hierFDR,setting)]
myTab3[hierFDR==T,table(phenotype,setting)]

save(myTab3, file = "../results/02_TWAS_SummaryStatistics_LIFEHeart_BMIadj_hierFDR.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

