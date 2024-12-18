#' ---
#' title: "MR - 2-sample - summary statistics"
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
#' - MR-IVW for each gene - hormone - setting combination
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "forostar"

source("../../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/07_MR_egger.RData")
load("../temp/eqtlgen_harmonized.RData")
load("../temp/UKB_IVs.RData")

#' # Add ENSG ID ####
#' ***
#' In the previous script, I forgot to add the ENSG ("Gene"). But I need this ID to ensure that I have unique IDs. 
#'  
MyGenes = eqtlgen_harmonized[,.N,by=c("Gene","GeneSymbol")]
dummy = MyGenes[,Gene]
dummy = rep(dummy,each=4)
MR_results[,ENSG := dummy]

#' # Hier FDR ####
#' ***
MR_results = MR_results[!is.na(pval_egger),]
myFDR = addHierarchFDR(pvalues = MR_results$pval_egger, categs = MR_results$exposure)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
MR_results[,P.Value.adj1 := myFDR$fdr_level1]
MR_results[,hierFDR := myFDR$hierarch_fdr5proz]
MR_results[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
MR_results[, table(hierFDR,exposure)]
MR_results[,min(P.Value.adj1),by=exposure]
MR_results[,min(pval_egger),by=exposure]
MR_results[,table(P.Value.adj1<0.05,exposure)]
MR_results[,table(pval_egger<0.05,exposure)]
MR_results[hierFDR == T,]

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

