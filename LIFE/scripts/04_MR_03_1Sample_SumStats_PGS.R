#' ---
#' title: "MR - 1-sample - summary statistics"
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
#' Plan: 
#' 
#' - TWAS of PGS on GE
#' - linear regression of PGS on hormones
#' - simple Wald ratio per GE 
#' - hierarchical FDR over Wald ratio p-values
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
#' Phenotypes 
#' 
load(paste0(path_LIFEprepped,"phenotypes/Adult_QC.RData"))
head(myTab)

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

myTab = myTab[eQTL==T,]
myTab[,table(group, is.na(daysLastMenst))]

#' Get GE ID
GE_samples = data.table(read_excel(path_LIFEAdult_GE_doku))
matched = match(myTab$ALIQUOT_GE,GE_samples$Aliquot)
myTab[,GE_ID := GE_samples[matched,GX_new_ID_v2]]

#' PGS 
#' 
load(paste0(path_LIFEprepped,"genetics/Adult_PGS_selection.RData"))
matched = match(myTab$ALIQUOT_genetics,myTab_Adult_PGS$Aliquot)
myTab = cbind(myTab,myTab_Adult_PGS[matched,-c(1,2)])
myTab[,CORT_comb := scale(CORT_comb)]
myTab[SEX == 1,TT_men := scale(TT_men)]
myTab[SEX == 1,E2_men := scale(E2_men)]
myTab[SEX == 2,TT_women := scale(TT_women)]

#' eSet 
#' 
loaded1 = load(path_LIFEAdult_GE)
loaded1
eset_A1

#' # Prep data for TWAS ####
#' ***
#' Add the relevant columns to the Expression Set
#' 
temp1 = pData(eset_A1)

toadd = c("TT_men","TT_women","E2_men","CORT_comb","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
          "activeSmoking","BMI","group","TimeBloodSample","daysLastMenst","OUremoved")
table(is.element(myTab$GE_ID,temp1$sampleID))
table(is.element(temp1$sampleID,myTab$GE_ID))

for(i in toadd) {
  pData(eset_A1)[i]    =  myTab[ match_hk(pData(eset_A1)$sampleID,  myTab$GE_ID), i, with = F]
}

eset_SH = eset_A1[,!is.na(eset_A1$CORT_comb)]
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
save(eset_SH, file = paste0(path_LIFEprepped,"transcriptomics/Adult_eSet_PGS.RData"))

#' # Run limma #### 
#' ***

#' ## Men ####
gx_assoc_men_filt = calculateLimma(todovars = c("TT_men","E2_men"), 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample",
                                            "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

gx_assoc_men_filt2= calculateLimma(todovars = c("TT_men","E2_men"), 
                                   conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI",
                                            "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                   eset = eset_SH[,eset_SH$group=="men"], 
                                   doscale = T) 

#' ## Women ####
gx_assoc_women_filt = calculateLimma(todovars = "TT_women", 
                                     conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group",
                                              "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                     eset = eset_SH[,eset_SH$group!="men"], 
                                     doscale = T) 

gx_assoc_women_filt2= calculateLimma(todovars = "TT_women", 
                                     conf = c("AGE","LYP","MOP","activeSmoking","TimeBloodSample","OUremoved","group","BMI",
                                              "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                     eset = eset_SH[,eset_SH$group!="men"], 
                                     doscale = T) 

#' ## Sex-combined ####
gx_assoc_comb_filt = calculateLimma(todovars = "CORT_comb", 
                                    conf = c("group","AGE","LYP","MOP","activeSmoking","TimeBloodSample",
                                             "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                    eset = eset_SH, 
                                    doscale = T) 

gx_assoc_comb_filt2= calculateLimma(todovars = "CORT_comb", 
                                    conf = c("group","AGE","LYP","MOP","activeSmoking","TimeBloodSample","BMI",
                                             "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                    eset = eset_SH, 
                                    doscale = T) 

#' # Extract Summary statistics
#' ***
#' I do not like lists - I want a data table!
#' 
ToDoList = data.table(hormones = c("TT","E2","TT","CORT"),
                      setting = c("men","men","women","comb"))
ToDoList[,PGS := paste(hormones,setting,sep="_")]
dumTab = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=4
  if(ToDoList[i,setting] == "women"){
    tab1 = myExtractionFunction(data=gx_assoc_women_filt,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_women_filt2,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
  }else if(ToDoList[i,setting] == "men"){
    tab1 = myExtractionFunction(data=gx_assoc_men_filt,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_men_filt2,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
  }else if(ToDoList[i,setting] == "comb"){
    tab1 = myExtractionFunction(data=gx_assoc_comb_filt,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
    tab2 = myExtractionFunction(data=gx_assoc_comb_filt2,SH = ToDoList[i,PGS],set = ToDoList[i,setting])
  }
  tab1[,model:="noBMIadj"]
  tab2[,model:="BMIadj"]
  tab = rbind(tab1,tab2)
  tab
}

Adult = rbindlist(dumTab)  
Adult[,study:="Adult"]
Adult[,uniqueID := paste(phenotype,model,PROBE_ID,sep="::")]
table(duplicated(Adult$uniqueID))

names(Adult)

Adult = Adult[,c(18,22,23,24,25,1,2,14,20,13,19,21,15:17,9:12,3:8)]
Adult

save(Adult, file = "../temp/04_MR_03_TWAS_SummaryStatistics_PGS.RData")

#' # Linear regression of PGS on hormones ####
#' ***
load("../results/02_PGS_10_PGSresults_selectionScaled.RData")
PGS.result = PGS.result[study=="Adult",]
stopifnot((PGS.result$score %in% Adult$phenotype))

#' # Wald estimator ####
#' ***
matched = match(Adult$phenotype,PGS.result$score)

Adult[,beta_GX := PGS.result[matched,beta]]
Adult[,SE_GX := PGS.result[matched,SE]]
Adult[,pval_GX := PGS.result[matched,pval]]

Adult[,beta_GY := beta]
Adult[,SE_GY := se]
Adult[,pval_GY := P.Value]

Adult[,beta_Wald := beta_GY/beta_GX]
Adult[,SE_Wald := sqrt(SE_GY^2/beta_GX^2 + beta_GY^2*SE_GX^2/beta_GX^4)]
Adult[,pval_Wald := 2*pnorm(-abs(beta_Wald/SE_Wald))]

table(Adult$pval_Wald<0.05,Adult$phenotype)

#' # Hierarchical FDR ####
#' ***
Adult[,trait := paste(phenotype,model,sep="_")]
Adult[,table(trait,pval_Wald<0.05)]

myFDR = addHierarchFDR(pvalues = Adult$pval_Wald, categs = Adult$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
Adult[,pval_Wald_adj1 := myFDR$fdr_level1]
Adult[,hierFDR := myFDR$hierarch_fdr5proz]
Adult[,pval_Wald_adj2 := myFDR$fdr_level1 * n / k]
Adult[, table(hierFDR,trait)]
Adult[,min(pval_Wald_adj1),by=c("phenotype","model")]


#' # Hierarchical FDR by gene set ####
#' ***
load("../results/03_TWAS_03_ProbeList2_Replicated_in_BMImodels.RData")
Adult2 = copy(Adult)
Adult2[,trait2 := paste(phenotype,PROBE_ID,sep="_")]
ProbeList2[phenotype=="TESTO",phenotype:="TT"]
ProbeList2[setting=="combined",setting:="comb"]
ProbeList2[,trait := paste(phenotype,setting,PROBE_ID,sep="_")]
Adult2 = Adult2[trait2 %in% ProbeList2$trait]
Adult2 = Adult2[model=="noBMIadj",]

myFDR = addHierarchFDR(pvalues = Adult2$pval_Wald, categs = Adult2$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
Adult2[,pval_Wald_adj1 := myFDR$fdr_level1]
Adult2[,hierFDR := myFDR$hierarch_fdr5proz]
Adult2[,pval_Wald_adj2 := myFDR$fdr_level1 * n / k]
Adult2[, table(hierFDR,trait)]
Adult2[,min(pval_Wald_adj1),by=c("phenotype","model")]

#' # Save results ####
#' ***
save(Adult,file = "../results/04_MR_03_1sample_SummaryStatistics_transcriptomeWide.RData")
save(Adult2,file = "../results/04_MR_03_1sample_SummaryStatistics_geneSet2.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

