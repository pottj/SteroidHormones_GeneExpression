#' ---
#' title: "MVMR - 2-sample - summary statistics"
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
#' - MVMR-IVW for each gene - setting combination
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
load("../temp/05_MVMR_eqtlgen_harmonized.RData")
load("../temp/05_MVMR_BMI_TT_IVs.RData")

#' # MR-IVW ####
#' ***
#' I need two loops: one over all genes (m=19873) and one over the hormones and sex-settings (n=4)
#' 
MyGenes = eqtlgen_harmonized[,.N,by=c("Gene","GeneSymbol")]
MyHormones = erg[,.N,by = c("setting")]

registerDoParallel(20)

dumTab3 = foreach(i = 1:dim(MyGenes)[1])%dopar%{
#dumTab3 = foreach(i = 1:10)%do%{
  #i=1
  myGene = copy(MyGenes)
  myGene = myGene[i,]
  #message("Working on gene ",myGene$GeneSymbol)
  
  outcome = copy(eqtlgen_harmonized)
  outcome = outcome[GeneSymbol == myGene$GeneSymbol,]
  
  dumTab4 = foreach(j = 1:dim(MyHormones)[1])%do%{
    #j=1
    myHormone = copy(MyHormones)
    myHormone = myHormone[j,]
    #message("   Working on hormone TT in ",myHormone$setting)
    
    exposure = copy(erg)
    exposure = exposure[setting == myHormone$setting,]
    
    outcome2 = copy(outcome)
    matched = match(exposure$rsID,outcome2$SNP)
    outcome2 = outcome2[matched]
    filt = is.na(outcome2$SNP)
    outcome2 = outcome2[!filt,]
    exposure = exposure[!filt,]
    stopifnot(outcome2$SNP == exposure$rsID)
    
    filt1 = grepl("BETA",names(exposure))
    data_beta = copy(exposure)
    data_beta = data_beta[,filt1,with=F]
    filt2 = grepl("SE",names(exposure))
    data_SE = copy(exposure)
    data_SE = data_SE[,filt2,with=F]
    types = gsub("BETA_","",names(exposure)[filt1])
    
    mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                          bxse = as.matrix(data_SE),
                          by = outcome2$beta, 
                          byse = outcome2$SE,
                          exposure = types,
                          outcome = myGene$GeneSymbol,
                          snps = exposure$rsID)
    
    res2 = mr_mvivw(mvmr_obj,nx = c(max(exposure$N_TT),max(exposure$N_BMI))) 
    n_SNPs_TT = dim(exposure[abs(Zscore_TT)>5,])[1]
    n_SNPs_BMI = dim(exposure[abs(Zscore_BMI)>5,])[1]
    
    res = data.table(setting = myHormone$setting,
                     exposure = c(res2@Exposure),
                     outcome = rep(res2@Outcome,length(types)),
                     NR_SNPs_total = rep(res2@SNPs,length(types)),
                     NR_SNPs_type = c(n_SNPs_TT,n_SNPs_BMI),
                     beta_IVW = c(res2@Estimate),
                     SE_IVW = c(res2@StdError),
                     pval_IVW = c(res2@Pvalue),
                     condFstat = c(res2@CondFstat),
                     HeteroStat = rep(res2@Heter.Stat[1],length(types)),
                     HeteroStat_pval = rep(res2@Heter.Stat[2],length(types)))
    
    res
    
  }
  MRIVW_results = rbindlist(dumTab4)
  MRIVW_results
}
MR_results = rbindlist(dumTab3)
MR_results[pval_IVW<0.05,]
save(MR_results,file = "../results/05_MVMR_02_IVW.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

