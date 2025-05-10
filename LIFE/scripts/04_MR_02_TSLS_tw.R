#' ---
#' title: "TSLS - gene set 1"
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

#' # Load data ####
#' ***
#' PGS 
#' 
load(paste0(path_LIFEprepped,"genetics/Adult_PGS_selection.RData"))
head(myTab_Adult_PGS)

#' eSet 
#' 
load(paste0(path_LIFEprepped,"transcriptomics/Adult_eSet.RData"))
dim(eset_SH)

#' Phenotypes 
#' 
load(paste0(path_LIFEprepped,"phenotypes/Adult_QC.RData"))
head(myTab)

#' Additional IDs
GE_samples = data.table(read_excel(path_LIFEAdult_GE_doku))
GE_samples = GE_samples[Aliquot %in% myTab$ALIQUOT_GE,]

#' # Match everything ####
#' ***
#' Match 1: alternative ID to phenotypes
#' 
matched = match(myTab$ALIQUOT_GE,GE_samples$Aliquot)
myTab[,GE_ID := GE_samples[matched,GX_new_ID_v2]]

#' Match 2: PGS to phenotypes
#' 
matched = match(myTab$ALIQUOT_genetics,myTab_Adult_PGS$Aliquot)
myTab = cbind(myTab,myTab_Adult_PGS[matched,-c(1,2)])

#' Filter 1: Samples with ID and hormone and PGS
#' 
myTab = myTab[!is.na(CORT) & !is.na(PC1) & !is.na(GE_ID),]
myTab[TSLS==F,]
myTab[, table(TSLS,is.na(CORT))]

#' Match 3: myTab and gene expression data
GEmat = eset_SH@assayData$exprs
dim(GEmat)

PhenoMat = eset_SH@phenoData@data
table(is.element(PhenoMat$sampleID,myTab$GE_ID))
filt = is.element(PhenoMat$sampleID,myTab$GE_ID)
PhenoMat = PhenoMat[filt,]
GEmat = GEmat[,filt]
myTab = myTab[GE_ID %in% PhenoMat$sampleID,]

matched = match(PhenoMat$sampleID,myTab$GE_ID)
myTab = myTab[matched,]
stopifnot(myTab$GE_ID == PhenoMat$sampleID)

#' # Loop ####
#' ***

counter = seq(1,20000,1000)
dumTab = foreach(i = 1:dim(GEmat)[1])%do%{
  #i=1
  #if(i %in% counter) message("Working on probe ",rownames(GEmat)[i],", number ",i," of ",dim(GEmat)[1])
  
  myGE = GEmat[i,]
  dat1 = copy(myTab)
  dat1[, GE := myGE]
  dat1[, CORT := log(CORT)]
  dat1[, TESTO := log(TESTO)]
  dat1[, E2 := log(E2)]
  
  filt_sex = dat1$GENDER==1
  
  ivmodel1 = ivreg(GE ~ CORT | CORT_comb, x=TRUE,data = dat1)
  ivmodel2 = ivreg(GE ~ TESTO | TT_men, x=TRUE,data = dat1[filt_sex, ])
  ivmodel3 = ivreg(GE ~ TESTO | TT_women, x=TRUE,data = dat1[!filt_sex, ])
  ivmodel4 = ivreg(GE ~ E2 | E2_men, x=TRUE,data = dat1[filt_sex, ])
  
  res = data.table(PROBE_ID = rep(rownames(GEmat)[i],4),
                   sampleSize = c(dim(dat1)[1],dim(dat1[filt_sex,])[1],dim(dat1[!filt_sex,])[1],dim(dat1[filt_sex,])[1]),
                   phenotype = c("CORT","TT","TT","E2"),
                   setting = c("comb","men","women","men"),
                   beta_2SLS = c(summary(ivmodel1)$coef[2,1],summary(ivmodel2)$coef[2,1],
                                 summary(ivmodel3)$coef[2,1],summary(ivmodel4)$coef[2,1]),
                   SE_2SLS = c(summary(ivmodel1)$coef[2,2],summary(ivmodel2)$coef[2,2],
                               summary(ivmodel3)$coef[2,2],summary(ivmodel4)$coef[2,2]),
                   tval_2SLS = c(summary(ivmodel1)$coef[2,3],summary(ivmodel2)$coef[2,3],
                                 summary(ivmodel3)$coef[2,3],summary(ivmodel4)$coef[2,3]),
                   pval_2SLS = c(summary(ivmodel1)$coef[2,4],summary(ivmodel2)$coef[2,4],
                                 summary(ivmodel3)$coef[2,4],summary(ivmodel4)$coef[2,4]))
  res
}
TSLS_Tab_Adult = rbindlist(dumTab)
TSLS_Tab_Adult[,min(pval_2SLS),by=c("phenotype","setting")]
save(TSLS_Tab_Adult,file = "../temp/04_MR_02_TSLS_transcriptomeWide.RData")

#' # Hierarchical FDR ####
#' ***
#' 
TSLS_Tab = copy(TSLS_Tab_Adult)
TSLS_Tab[,trait := paste(phenotype,setting,sep="_")]
#TSLS_Tab = TSLS_Tab[trait %in% c("TT_men","TT_women","CORT_comb"),]

myFDR = addHierarchFDR(pvalues = TSLS_Tab$pval_2SLS, categs = TSLS_Tab$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab[, table(hierFDR,trait)]
TSLS_Tab[,min(P.Value.adj1),by=c("phenotype","setting")]

#' # Save results ####
#' ***
save(TSLS_Tab,file = "../results/04_MR_02_TSLS_transcriptomeWide.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

