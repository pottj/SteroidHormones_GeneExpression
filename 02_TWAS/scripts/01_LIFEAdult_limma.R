#' ---
#' title: "TWAS - LIFE-Adult"
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
server = "angmar"

source("../../SourceFile.R")
.libPaths()

source(func_TWAS)
source(func_TWASExtract)

#' # Load data ####
#' ***
#' Load the phenotype data (hormones and covariates) and the Expression Set (GE)
#' 
load(paste0(path_LIFEprepped,"phenotypes/Adult_QC.RData"))
myTab

loaded1 = load(path_LIFEAdult_GE)
loaded1
eset_A1

GE_samples = data.table(read_excel(path_LIFEAdult_GE_doku))
GE_samples = GE_samples[Aliquot %in% myTab$ALIQUOT_GE,]
matched = match(myTab$ALIQUOT_GE,GE_samples$Aliquot)
myTab[,GE_ID := GE_samples[matched,GX_new_ID_v2]]

#' # Transform data ####
#' ***
#' I want to log-transform the hormone data and BMI, make smoking binary (active smoking only), and change some column names
#' 
names(myTab)
myTab[,CORT := log(CORT)]
myTab[,TESTO := log(TESTO)]
myTab[,E2 := log(E2)]
myTab[,D074_BMI := log(D074_BMI)]

myTab[,table(D141_smokeStatus)]
myTab[D141_smokeStatus==1,D141_smokeStatus:=0]
myTab[D141_smokeStatus==2,D141_smokeStatus:=1]

setnames(myTab,"GENDER","SEX")
setnames(myTab,"D126_time","TimeBloodSample")
setnames(myTab,"D126_fasting","FastingHours")
setnames(myTab,"D074_BMI","BMI")
setnames(myTab,"D141_smokeStatus","activeSmoking")
setnames(myTab,"D133_daysLastMenst","daysLastMenst")

myTab = myTab[TWAS==T,]
myTab[,table(group, is.na(daysLastMenst))]

#' # Merge data ####
#' ***
#' Add the relevant columns to the Expression Set
#' 
temp1 = pData(eset_A1)

toadd = c("CORT","TESTO","E2" ,
          "activeSmoking","BMI","group","TimeBloodSample","daysLastMenst","OUremoved")
table(is.element(myTab$GE_ID,temp1$sampleID))
table(is.element(temp1$sampleID,myTab$GE_ID))

for(i in toadd) {
  pData(eset_A1)[i]    =  myTab[ match_hk(pData(eset_A1)$sampleID,  myTab$GE_ID), i, with = F]
}

#' # Filter data ####
#' ***
#' Reduce Expression Set to samples with hormone data 
eset_SH = eset_A1[,!is.na(eset_A1$CORT)]
dim(eset_SH)
temp2 = pData(eset_SH)
table(is.na(temp2$CORT),temp2$group)

#' Select high quality probes only
temp3 = fData(eset_SH)
filt_features<-temp3$probe_QCok==T & temp3$unique_mapper==T
table(filt_features)

goodfeatures<-temp3$PROBE_ID[filt_features]
eset_SH = eset_SH[goodfeatures,]
dim(eset_SH)

#' save filtered eSet
#' 
save(eset_SH, file = paste0(path_LIFEprepped,"transcriptomics/Adult_eSet.RData"))

#' # Run limma #### 
#' ***
SHs_sex = c("TESTO","E2")
SHs_comb = c("CORT")

#' ## Men ####
gx_assoc_men_filt = calculateLimma(todovars = SHs_sex, 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

gx_assoc_men_filt2= calculateLimma(todovars = SHs_sex, 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

#' ## Women ####
gx_assoc_women_filt = calculateLimma(todovars = SHs_sex, 
                                        conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group"), 
                                        eset = eset_SH[,eset_SH$group!="men"], 
                                        doscale = T) 

gx_assoc_women_filt2= calculateLimma(todovars = SHs_sex, 
                                        conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group","BMI"), 
                                        eset = eset_SH[,eset_SH$group!="men"], 
                                        doscale = T) 

#' ## Sex-combined ####
gx_assoc_comb_filt = calculateLimma(todovars = SHs_comb, 
                                     conf = c("group","AGE","LYP","MOP","activeSmoking","TimeBloodSample"), 
                                     eset = eset_SH, 
                                     doscale = T) 

gx_assoc_comb_filt2= calculateLimma(todovars = SHs_comb, 
                                     conf = c("group","AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI"), 
                                     eset = eset_SH, 
                                     doscale = T) 

#' # Save lists ####
#' ***
save(gx_assoc_men_filt, file="../temp/01_TWAS_Adult_men.RData")
save(gx_assoc_men_filt2, file="../temp/01_TWAS_Adult_men_BMIadj.RData")

save(gx_assoc_women_filt, file="../temp/01_TWAS_Adult_women.RData")
save(gx_assoc_women_filt2, file="../temp/01_TWAS_Adult_women_BMIadj.RData")

save(gx_assoc_comb_filt, file="../temp/01_TWAS_Adult_comb.RData")
save(gx_assoc_comb_filt2, file="../temp/01_TWAS_Adult_comb_BMIadj.RData")

#' # Extract Summary statistics
#' ***
#' I do not like lists - I want a data table!
#' 
ToDoList = data.table(hormones = c(SHs_sex,SHs_sex,SHs_comb),
                      setting = c("men","men","women","women","combined"))

dumTab = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  if(ToDoList[i,setting] == "women"){
    tab1 = myExtractionFunction(data=gx_assoc_women_filt,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_women_filt2,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
  }else if(ToDoList[i,setting] == "men"){
    tab1 = myExtractionFunction(data=gx_assoc_men_filt,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_men_filt2,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
  }else if(ToDoList[i,setting] == "combined"){
    tab1 = myExtractionFunction(data=gx_assoc_comb_filt,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_comb_filt2,SH = ToDoList[i,hormones],set = ToDoList[i,setting])
  }
  tab1[,model:="noBMIadj"]
  tab2[,model:="BMIadj"]
  tab = rbind(tab1,tab2)
  tab
}

Adult = rbindlist(dumTab)  
Adult[,study:="Adult"]
Adult[,uniqueID := paste(phenotype,setting,model,PROBE_ID,sep="::")]
table(duplicated(Adult$uniqueID))

names(Adult)

Adult = Adult[,c(18,22,23,24,25,1,2,14,20,13,19,21,15:17,9:12,3:8)]
Adult

save(Adult, file = "../results/01_TWAS_SummaryStatistics_LIFEAdult.RData")

#' # Hierarchical FDR ####
#' ***
myTab = copy(Adult)
myTab[,table(setting,phenotype,model)]
myTab[,trait := paste(phenotype,setting,model,sep="_")]
myTab[,table(trait)]

myFDR = addHierarchFDR(pvalues = myTab$P.Value, categs = myTab$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
myTab[,P.Value.adj1 := myFDR$fdr_level1]
myTab[,hierFDR := myFDR$hierarch_fdr5proz]
myTab[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
myTab[, table(hierFDR,trait)]

save(myTab, file = "../results/01_TWAS_SummaryStatistics_LIFEAdult_hierFDR.RData")

#' # Save per trait ####
#' ***
#' I plan to upload the TWAS summary statistics per trait. Hence, I want a text file per hormone - setting - model
#' 
outdir = "../results/01_Adult_SummaryStatistics/"

myTraits = unique(myTab$trait)

for(i in 1:length(myTraits)){
  #i=1
  myTab2 = copy(myTab)
  myTab2 = myTab2[trait == myTraits[i],]
  myTab2 = myTab2[,c(1:3,6:14,27,29,28,16:25)]
  setorder(myTab2,ilmn_chr_hg19,ilmn_start_hg19)
  
  outfn = paste0(outdir,myTraits[i],".txt")
  fwrite(myTab2, file = outfn,quote = F,sep="\t")
  gzip(outfn,destname = paste0(outfn, ".gz"))
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

