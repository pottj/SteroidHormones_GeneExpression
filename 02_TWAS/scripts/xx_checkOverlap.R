#' ---
#' title: "Check overlap betweeen LIFE"
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

#' # Load data ####
#' ***
loaded1 = load("../results/01_TWAS_SummaryStatistics_LIFEAdult_hierFDR.RData")
myTabA1 = get(loaded1)
loaded2 = load("../results/02_TWAS_SummaryStatistics_LIFEHeart_hierFDR.RData")
myTabB3 = get(loaded2)

#' # Filter ####
#' ***
SharedProbes = myTabA1[hierFDR == T & PROBE_ID %in% myTabB3[hierFDR==T,PROBE_ID], unique(PROBE_ID)]

myTabA1 = myTabA1[PROBE_ID %in% SharedProbes,]
myTabB3 = myTabB3[PROBE_ID %in% SharedProbes,]

myTabA1[hierFDR==T & uniqueID %in% myTabB3[hierFDR==T,uniqueID],table(phenotype,setting)]

#' # Cortisol ####
#' ***
myTabA1 = myTabA1[hierFDR==T & uniqueID %in% myTabB3[hierFDR==T,uniqueID],]
myTabB3 = myTabB3[hierFDR==T & uniqueID %in% myTabA1[hierFDR==T,uniqueID],]

myTabA1_cort = copy(myTabA1)
myTabA1_cort = myTabA1_cort[phenotype=="CORT"]

myTabB3_cort = copy(myTabB3)
myTabB3_cort = myTabB3_cort[phenotype=="CORT"]

matched = match(myTabA1_cort$uniqueID,myTabB3_cort$uniqueID)
myTabB3_cort = myTabB3_cort[matched,]

table(myTabB3_cort$setting == myTabA1_cort$setting)

plot(myTabA1_cort$beta,myTabB3_cort$beta)
abline(0,1)

test = myTabA1_cort[,.N,by = PROBE_ID]
setorder(myTabA1_cort,P.Value.adj1)
myTabA1_cort = myTabA1_cort[!duplicated(PROBE_ID),]
myTabA1_cort[,table(setting)]
