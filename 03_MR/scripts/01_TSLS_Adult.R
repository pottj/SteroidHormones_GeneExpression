#' ---
#' title: "TSLS"
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

source("../../SourceFile_angmar.R")
.libPaths()

#' # Load data ####
#' ***
#' Filtered eset and phenotype file with PGSs
#' 
load(paste0(path_LIFEprepped,"05_LIFEAdult_filtered_final_PGS.RData"))
myTab = myTab[TSLS == T,]

load(paste0(path_LIFEprepped,"eSet_Adult_TWAS.RData"))
dim(eset_SH)

GE_samples = data.table(read_excel(path_LIFEAdult_GE_doku))
GE_samples = GE_samples[Aliquot %in% myTab$ALIQUOT_GE,]
matched = match(myTab$ALIQUOT_GE,GE_samples$Aliquot)
myTab[,GE_ID := GE_samples[matched,GX_new_ID_v2]]

GEmat = eset_SH@assayData$exprs
dim(GEmat)

PhenoMat = eset_SH@phenoData@data
table(is.element(PhenoMat$sampleID,myTab$GE_ID))
filt = is.element(PhenoMat$sampleID,myTab$GE_ID)
PhenoMat = PhenoMat[filt,]
GEmat = GEmat[,filt]

matched = match(PhenoMat$sampleID,myTab$GE_ID)
myTab = myTab[matched,]
stopifnot(myTab$GE_ID == PhenoMat$sampleID)

#' # Loop ####
#' ***
load("../../02_TWAS/results/01_TWAS_SummaryStatistics_LIFEAdult_hierFDR.RData")
myTab2 = myTab2[hierFDR==T,]
filt = is.element(rownames(GEmat),myTab2$PROBE_ID)
table(filt)
GEmat = GEmat[filt,]

counter = seq(1,10000,100)
dumTab = foreach(i = 1:dim(GEmat)[1])%do%{
  #i=1
  if(i %in% counter) message("Working on probe ",i)
  
  myGE = GEmat[i,]
  dat1 = copy(PhenoMat)
  dat1$GE = myGE
  dat1$PGS_TT = myTab$PGS_TT_same
  dat1$PGS_E2 = myTab$PGS_E2_same
  dat1$PGS_CORT = myTab$PGS_CORT_Combined_best
  
  filt = dat1$SEX==1
  ivmodel.TT.men = ivreg(GE~TESTO|PGS_TT, x=TRUE,data = dat1[filt,])
  ivmodel.TT.women = ivreg(GE~TESTO|PGS_TT, x=TRUE,data = dat1[!filt,])
  ivmodel.E2.men = ivreg(GE~E2|PGS_E2, x=TRUE,data = dat1[filt,])
  ivmodel.E2.women = ivreg(GE~E2|PGS_E2, x=TRUE,data = dat1[!filt,])
  ivmodel.CORT = ivreg(GE~CORT|PGS_CORT, x=TRUE,data = dat1)

  res = data.table(PROBE_ID = rownames(GEmat)[i],
                   sampleSize = c(1296, 1005, 1296, 1005, 2301),
                   phenotype = c("TT","TT","E2","E2","CORT"),
                   sex = c("men","women","men","women","combined"),
                   beta_2SLS = c(summary(ivmodel.TT.men)$coef[2,1],summary(ivmodel.TT.women)$coef[2,1],
                                 summary(ivmodel.E2.men)$coef[2,1],summary(ivmodel.E2.women)$coef[2,1],
                                 summary(ivmodel.CORT)$coef[2,1]),
                   SE_2SLS = c(summary(ivmodel.TT.men)$coef[2,2],summary(ivmodel.TT.women)$coef[2,2],
                               summary(ivmodel.E2.men)$coef[2,2],summary(ivmodel.E2.women)$coef[2,2],
                               summary(ivmodel.CORT)$coef[2,2]),
                   tval_2SLS = c(summary(ivmodel.TT.men)$coef[2,3],summary(ivmodel.TT.women)$coef[2,3],
                                 summary(ivmodel.E2.men)$coef[2,3],summary(ivmodel.E2.women)$coef[2,3],
                                 summary(ivmodel.CORT)$coef[2,3]),
                   pval_2SLS = c(summary(ivmodel.TT.men)$coef[2,4],summary(ivmodel.TT.women)$coef[2,4],
                                 summary(ivmodel.E2.men)$coef[2,4],summary(ivmodel.E2.women)$coef[2,4],
                                 summary(ivmodel.CORT)$coef[2,4]))
  res
}
TSLS_Tab_Adult = rbindlist(dumTab)
TSLS_Tab_Adult[,min(pval_2SLS),by=c("phenotype","sex")]

test = copy(TSLS_Tab_Adult)
test = test[pval_2SLS<0.05,]

#' # Save results ####
#' ***
save(TSLS_Tab_Adult,file = "../results/01_TSLS_Adult.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

