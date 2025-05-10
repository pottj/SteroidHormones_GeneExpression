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

#' Candidate genes
load("../results/03_TWAS_03_ProbeList5_Adult_any.RData")

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

#' Filter 2: select relevant genes
#' 
ProbeList = ProbeList5[!grepl("E2_women",trait),]
ProbeList[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
table(duplicated(ProbeList$trait2))
ProbeList = ProbeList[!duplicated(trait2),]
length(unique(ProbeList$PROBE_ID))

filt = is.element(rownames(GEmat),ProbeList$PROBE_ID)
table(filt)
GEmat = GEmat[filt,]
dim(GEmat)

#' # Loop ####
#' ***

counter = seq(1,10000,100)
dumTab = foreach(i = 1:dim(ProbeList)[1])%do%{
  #i=1
  myProbe = ProbeList[i,]
  #if(i %in% counter) message("Working on probe ",myProbe$PROBE_ID,", number ",i," of ",dim(ProbeList)[1])
  
  x = grep(myProbe$PROBE_ID,rownames(GEmat))
  myGE = GEmat[x,]
  dat1 = copy(myTab)
  dat1[, GE := myGE]
  dat1[, hormone := get(myProbe$phenotype)]
  dat1[, hormone := log(hormone)]
  
  if(myProbe$setting=="men"){
    dat1 = dat1[group=="men"]
  }else if(myProbe$setting=="women"){
    dat1 = dat1[group!="men"]
  }
  
  if(myProbe$phenotype=="CORT"){
    dat1[,PGS := CORT_comb]
  }else if(myProbe$phenotype=="TESTO"){
    if(myProbe$setting=="men"){
      dat1[,PGS := TT_men]
    }else if(myProbe$setting=="women"){
      dat1[,PGS := TT_women]
    }
  }else if(myProbe$phenotype=="E2"){
    if(myProbe$setting=="men"){
      dat1[,PGS := E2_men]
    }else if(myProbe$setting=="women"){
      dat1[,PGS := E2_women]
    }
  }
  ivmodel = ivreg(GE~hormone|PGS, x=TRUE,data = dat1)
    
  summary(ivmodel)
  res = data.table(PROBE_ID = myProbe$PROBE_ID,
                   sampleSize = dim(dat1)[1],
                   phenotype = myProbe$phenotype,
                   setting = myProbe$setting,
                   beta_2SLS = c(summary(ivmodel)$coef[2,1]),
                   SE_2SLS = c(summary(ivmodel)$coef[2,2]),
                   tval_2SLS = c(summary(ivmodel)$coef[2,3]),
                   pval_2SLS = c(summary(ivmodel)$coef[2,4]))
  res
}
TSLS_Tab_Adult = rbindlist(dumTab)
TSLS_Tab_Adult[,min(pval_2SLS),by=c("phenotype","setting")]

#' # Hierarchical FDR ####
#' ***
#' I decide to exclude the E2 in men and TT in women, as they are unlikely to do anything
#' 
#' ## Gene set 5 ####
TSLS_Tab_5 = copy(TSLS_Tab_Adult)
TSLS_Tab_5[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab_5 = TSLS_Tab_5[trait %in% c("TESTO_men","CORT_combined"),]

myFDR = addHierarchFDR(pvalues = TSLS_Tab_5$pval_2SLS, categs = TSLS_Tab_5$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab_5[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab_5[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab_5[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab_5[, table(hierFDR,trait)]
TSLS_Tab_5[,min(P.Value.adj1),by=c("phenotype","setting")]

#' ## Gene set 4 ####
load("../results/03_TWAS_03_ProbeList4_Adult_BMIadj.RData")
TSLS_Tab_4 = copy(TSLS_Tab_Adult)
TSLS_Tab_4[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab_4 = TSLS_Tab_4[trait %in% c("TESTO_men","CORT_combined"),]
TSLS_Tab_4[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
ProbeList4[,trait := paste(phenotype,setting,PROBE_ID,sep="_")]
TSLS_Tab_4 = TSLS_Tab_4[trait2 %in% ProbeList4$trait]

myFDR = addHierarchFDR(pvalues = TSLS_Tab_4$pval_2SLS, categs = TSLS_Tab_4$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab_4[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab_4[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab_4[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab_4[, table(hierFDR,trait)]
TSLS_Tab_4[,min(P.Value.adj1),by=c("phenotype","setting")]

#' ## Gene set 3 ####
load("../results/03_TWAS_03_ProbeList3_Adult_noBMIadj.RData")
TSLS_Tab_3 = copy(TSLS_Tab_Adult)
TSLS_Tab_3[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab_3 = TSLS_Tab_3[trait %in% c("TESTO_men","CORT_combined"),]
TSLS_Tab_3[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
ProbeList3[,trait := paste(phenotype,setting,PROBE_ID,sep="_")]
TSLS_Tab_3 = TSLS_Tab_3[trait2 %in% ProbeList3$trait]

myFDR = addHierarchFDR(pvalues = TSLS_Tab_3$pval_2SLS, categs = TSLS_Tab_3$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab_3[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab_3[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab_3[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab_3[, table(hierFDR,trait)]
TSLS_Tab_3[,min(P.Value.adj1),by=c("phenotype","setting")]

#' ## Gene set 2 ####
load("../results/03_TWAS_03_ProbeList2_Replicated_in_BMImodels.RData")
TSLS_Tab_2 = copy(TSLS_Tab_Adult)
TSLS_Tab_2[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab_2 = TSLS_Tab_2[trait %in% c("TESTO_men","CORT_combined"),]
TSLS_Tab_2[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
ProbeList2[,trait := paste(phenotype,setting,PROBE_ID,sep="_")]
TSLS_Tab_2 = TSLS_Tab_2[trait2 %in% ProbeList2$trait]

myFDR = addHierarchFDR(pvalues = TSLS_Tab_2$pval_2SLS, categs = TSLS_Tab_2$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab_2[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab_2[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab_2[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab_2[, table(hierFDR,trait)]
TSLS_Tab_2[,min(P.Value.adj1),by=c("phenotype","setting")]

#' ## Gene set 2 ####
load("../results/03_TWAS/03_ProbeList1_Replicated_in_Heart.RData")
TSLS_Tab_1 = copy(TSLS_Tab_Adult)
TSLS_Tab_1[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab_1 = TSLS_Tab_1[trait %in% c("TESTO_men","CORT_combined"),]
TSLS_Tab_1[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
ProbeList1[,trait := paste(phenotype,setting,PROBE_ID,sep="_")]
TSLS_Tab_1 = TSLS_Tab_1[trait2 %in% ProbeList1$trait]

myFDR = addHierarchFDR(pvalues = TSLS_Tab_1$pval_2SLS, categs = TSLS_Tab_1$trait)
k = length(unique(myFDR[hierarch_fdr5proz==T,category]))
n = length(unique(myFDR[,category]))
table(myFDR$hierarch_fdr5proz)
TSLS_Tab_1[,P.Value.adj1 := myFDR$fdr_level1]
TSLS_Tab_1[,hierFDR := myFDR$hierarch_fdr5proz]
TSLS_Tab_1[,P.Value.adj2 := myFDR$fdr_level1 * n / k]
TSLS_Tab_1[, table(hierFDR,trait)]
TSLS_Tab_1[,min(P.Value.adj1),by=c("phenotype","setting")]

#' # Save results ####
#' ***
#' I want to add the adjusted p-values per gene set, and then save the data table.
TSLS_Tab = copy(TSLS_Tab_Adult)
TSLS_Tab[,trait := paste(phenotype,setting,sep="_")]
TSLS_Tab = TSLS_Tab[trait %in% c("TESTO_men","CORT_combined"),]
TSLS_Tab[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]
TSLS_Tab_5[,trait2 := paste(phenotype,setting,PROBE_ID,sep="_")]

matched = match(TSLS_Tab$trait2,TSLS_Tab_1$trait2)
TSLS_Tab[,P.Value.adj1_GeneSet1 := TSLS_Tab_1[matched,P.Value.adj1]]

matched = match(TSLS_Tab$trait2,TSLS_Tab_2$trait2)
TSLS_Tab[,P.Value.adj1_GeneSet2 := TSLS_Tab_2[matched,P.Value.adj1]]

matched = match(TSLS_Tab$trait2,TSLS_Tab_3$trait2)
TSLS_Tab[,P.Value.adj1_GeneSet3 := TSLS_Tab_3[matched,P.Value.adj1]]

matched = match(TSLS_Tab$trait2,TSLS_Tab_4$trait2)
TSLS_Tab[,P.Value.adj1_GeneSet4 := TSLS_Tab_4[matched,P.Value.adj1]]

matched = match(TSLS_Tab$trait2,TSLS_Tab_5$trait2)
TSLS_Tab[,P.Value.adj1_GeneSet5 := TSLS_Tab_5[matched,P.Value.adj1]]

save(TSLS_Tab,file = "../results/04_MR_01_TSLS_Adult.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

