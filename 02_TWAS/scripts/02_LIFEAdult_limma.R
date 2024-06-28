#' ---
#' title: "Run limma in LIFE-Adult"
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
load(paste0(path_LIFEprepped,"01_LIFEAdult_filtered_final.RData"))
myTab

data = list.files(path = path_LIFEAdult, pattern = ".xlsx")
D00403 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00403",data)])))

loaded1 = load(path_LIFEAdult_GE)
loaded1

eset_A1

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
D00403 = D00403[D00403_ALIQUOT %in% myTab$ALIQUOT_GE,]

#' # Merge data ####
#' ***
#' Add the relevant columns to the Expression Set
#' 
temp1 = pData(eset_A1)

toadd = c("CORT","TESTO","E2" ,"activeSmoking","BMI","group","TimeBloodSample","daysLastMenst")
table(is.element(temp1$Aliquot,myTab$ALIQUOT_SH))

for(i in toadd) {
  pData(eset_A1)[i]    =  myTab[ match_hk(pData(eset_A1)$Aliquot,  myTab$ALIQUOT_SH), i, with = F]
}

#' # Filter data ####
#' ***
#' Reduce Expression Set to samples with hormone data & high quality probes
#' Sample-Filter: keep best sample in case of duplicates
eset_SH = eset_A1[,eset_A1$BEST_SAMPLE==T]
dim(eset_SH)
eset_SH = eset_A1[,!is.na(eset_A1$CORT)]
dim(eset_SH)

#' 
#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

