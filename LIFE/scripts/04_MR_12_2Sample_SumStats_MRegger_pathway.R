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
load("../temp/04_MR_eqtlgen_harmonized_pathway.RData")
load("../temp/04_MR_UKB_IVs_pathway.RData")

#' # MR-IVW ####
#' ***
#' I need two loops: one over all genes (m=19873) and one over the hormones and sex-settings (n=4)
#' 
MyGenes = eqtlgen_harmonized[,.N,by=c("Gene","GeneSymbol")]
MyHormones = erg[,.N,by = c("phenotype","setting")]

dumTab3 = foreach(i = 1:dim(MyGenes)[1])%do%{
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
    #message("   Working on hormone ",myHormone$phenotype," in ",myHormone$setting)
    
    exposure = copy(erg)
    exposure = exposure[phenotype == myHormone$phenotype & setting == myHormone$setting,]
    
    outcome2 = copy(outcome)
    matched = match(exposure$rsID,outcome2$SNP)
    outcome2 = outcome2[matched]
    filt = is.na(outcome2$SNP)
    outcome2 = outcome2[!filt,]
    exposure = exposure[!filt,]
    stopifnot(outcome2$SNP == exposure$rsID)
    
    if(dim(outcome2)[1]>2){
      mr_obj = mr_input(bx = as.numeric(exposure$BETA),
                        bxse = as.numeric(exposure$SE),
                        by = as.numeric(outcome2$beta),
                        byse = as.numeric(outcome2$SE),
                        exposure = paste(myHormone$phenotype,myHormone$setting,sep="_"),
                        outcome = myGene$GeneSymbol,
                        snps = exposure$rsID)
      res2 = mr_egger(mr_obj)
      res3 = mr_ivw(mr_obj)
      #mr_plot(mr_obj,line = "egger", orientate = T)
      res = data.table(exposure = res2@Exposure,
                       outcome = res2@Outcome,
                       NR_SNPs_total = res2@SNPs,
                       FStat = res3@Fstat,
                       beta_egger = res2@Estimate,
                       SE_egger = res2@StdError.Est,
                       pval_egger = res2@Pvalue.Est,
                       beta_int = res2@Intercept,
                       SE_int = res2@StdError.Int,
                       pval_int = res2@Pvalue.Int,
                       HeteroStat = res2@Heter.Stat[1],
                       HeteroStat_pval = res2@Heter.Stat[2])
      
    }else{
      mr_obj = mr_input(bx = as.numeric(exposure$BETA),
                        bxse = as.numeric(exposure$SE),
                        by = as.numeric(outcome2$beta),
                        byse = as.numeric(outcome2$SE),
                        exposure = paste(myHormone$phenotype,myHormone$setting,sep="_"),
                        outcome = myGene$GeneSymbol,
                        snps = exposure$rsID)
      res3 = mr_ivw(mr_obj)
      res = data.table(exposure = res3@Exposure,
                       outcome = res3@Outcome,
                       NR_SNPs_total = res3@SNPs,
                       comment = "not enough SNPs for MR-Egger")
    }
    
    res
    
  }
  MRIVW_results = rbindlist(dumTab4,fill=T,use.names = T)
  MRIVW_results
}
MR_results = rbindlist(dumTab3,fill=T,use.names = T)
MR_results[pval_egger<0.05,]
save(MR_results,file = "../results/04_MR_12_MR_egger_pathway.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

