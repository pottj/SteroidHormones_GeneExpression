#' ---
#' title: "Run 2SLS analyses"
#' subtitle: "MR - Testosterone in men on protein levels (Olink)"
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
#' The screening analysis is simply a 2SLS approach: 
#' 
#' protein ~ testosterone | 9 SNPs within/near genes coding for enzyms of the steroid hormone pathway
#' 
#' I want to restrict the analysis to proteins for which I also have gene expression data (eQTLGen phase I). I think it will be interesting to compare these findings later. 
#' 
#' I will add to my gene list the information if they are regulated by an androgen response element (ARE) in prostate cell lines (Wilson et al., 2016), and if they were associated in my TWAS in LIFE-Adult (GE in whole blood)
#' 
#' If there are any significant results, I will attempt a 2-sample MR using Sun et al. for the proteomics data, and Ruth et al. for the testosterone data (higher power because higher sample sizes?)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

library(data.table)
setDTthreads(1)
library(pgenlibr)
library(foreach)
library(doParallel)
library(ivreg)

#' # Parameter settings ####
#' ***
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"
UKB_proteomics = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics/"

#' # Load data ####
#' ***
#' Load phenotype file 
load(paste0(data_QC,"/01_Prep_01_UKB.RData"))

#' Load protein info 
load(paste0(data_QC,"/01_Prep_03_proteinIDs.RData"))

#' Load genotypes
data_QC2 = gsub("~","../../",data_QC)
pvar1 = NewPvar(paste0(data_QC2, '01_Prep_02_UKB_Testo_merged.pvar'))
pgen = NewPgen(paste0(data_QC2,'01_Prep_02_UKB_Testo_merged.pgen'), pvar=pvar1)
pvar = fread(paste0(data_QC2,'01_Prep_02_UKB_Testo_merged.pvar'))
psam = fread(paste0(data_QC2,'01_Prep_02_UKB_Testo_merged.psam'))

setnames(psam, "#FID","FID")
setnames(pvar, "#CHROM","CHR")

myNRs = 1:dim(pvar)[1]
geno_mat = ReadList(pgen, variant_subset = myNRs , meanimpute=F)
dim(geno_mat)
colnames(geno_mat) = pvar[,ID]
rownames(geno_mat) = psam[,IID]

#' Load Olink data
olink = fread(paste0(UKB_proteomics,"olink_data.txt.gz"))

#' # Filter Olink data ####
#' ***
#' Only samples I have in myTab & only probes I have in myGenes
#' 
olink = olink[eid %in% myTab$ID,]
olink = olink[ins_index == 0,]
olink = olink[protein_id %in% myGenes$UKB_coding,]

#' # Merge genotypes to phenotypes ####
#' ***
table(psam$SEX)
table(is.element(psam$IID,myTab$ID))
matched1 = match(myTab$ID,psam$IID)
pvar
myTab = cbind(myTab,geno_mat[matched1,])
myTab[,geneticSex := psam[matched1,SEX]]

#' # Loop per protein ####
#' ***
counter = seq(1,dim(myGenes)[1],50)
registerDoParallel(10)

dumTab1 = foreach(i = 1:dim(myGenes)[1])%dopar%{
  #i=1
  myRow = myGenes[i,]
  if(i %in% counter) message("Working on gene ",myRow$GeneSymbol,
                             " (number ",i," of ",dim(myGenes)[1],")")
  
  # filter Olink data
  data1 = copy(olink)
  data1 = data1[protein_id == myRow$UKB_coding,]
  
  # add to main data table
  matched2 = match(myTab$ID,data1$eid)
  myTab[,protein := data1[matched2,result]]
  
  # 2stage least square regression
  ivmodel = ivreg(protein ~ TESTO | rs559555 + rs11563251 + rs7694379 + rs4646450 + rs1832007 + rs28929474 + rs2414095 + rs12150660, x=TRUE,data = myTab)
  
  summary(ivmodel)
  res = data.table(exposure = "TESTO",
                   outcome = myRow$GeneSymbol,
                   sampleSize = dim(myTab[!is.na(protein)])[1],
                   beta_2SLS = c(summary(ivmodel)$coef[2,1]),
                   SE_2SLS = c(summary(ivmodel)$coef[2,2]),
                   tval_2SLS = c(summary(ivmodel)$coef[2,3]),
                   pval_2SLS = c(summary(ivmodel)$coef[2,4]))
  res
}
TSLS_Tab = rbindlist(dumTab1)
TSLS_Tab[,table(pval_2SLS<0.05)]

save(TSLS_Tab,file = "../results/02_2SLS_results.RData")

myTab[,protein := NULL]
save(myTab,file = paste0(data_QC,"/02_2SLS_input1.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
