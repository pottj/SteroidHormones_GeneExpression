#############################
# this is a template source file
# please change all parameters accordingly
#############################
# make files readable for all
Sys.umask(mode = "0002") 

#############################
# R library and R packages
#############################
# used R versions: R >= 4.2.x

if(server=="angmar"){
  .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/angmar/")
  suppressPackageStartupMessages(library(ivreg))
  suppressPackageStartupMessages(library(pgenlibr))
}else if(server=="forostar"){
  .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/forostar/")
  suppressPackageStartupMessages(library(MendelianRandomization))
}

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(toolboxH))

basicpath = "/net/ifs1/san_projekte/projekte/genstat/"
projectpath = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/2406_lifea1b3_MR_SH_GE/"

#############################
# Other Tools 
#############################
# PLINK2 path
path_plink2 = paste0(basicpath,"/07_programme/plink2.0/20210203/unix_64/plink2")
path_plink1.9 = paste0(basicpath,"/07_programme/plink1.9/vs210606_stable_beta_6_24/unix64/plink")

#############################
# Other helper functions 
#############################
# some of these functions are Holgers - I cannot share without his consent
func_TWAS = paste0(basicpath,"07_programme/rtools/1807_gx_tools/AllgFunktionenExpressionspipeline_220927.R") 

# these are my functions
func_TWASExtract = "../../helperfunctions/myExtractionFunction.R"


#############################
# path to LIFE data 
#############################
# all individual level data is not part of the repository!
#
# phenotype data
path_LIFEAdult = paste0(projectpath,"_data/pv785/data_adult/")
path_LIFEHeart = paste0(projectpath,"_data/pv785/data_heart/")
path_LIFEHeart_old = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1910_lifea1b3_sasha/_data/2020-03-30_PV505/heart/"

# gene expression data
path_LIFEAdult_GE = paste0(basicpath,"/02_qcdaten/1905_LifeA1_gx/ExpressionSet_Prepro_v3_LIFE_Adult_probesUNfiltered_IPA-2020-09-24.RData")
path_LIFEAdult_GE_doku = paste0(basicpath,"/02_qcdaten/1905_LifeA1_gx/Doku_Tab_LIFE_Adult.xlsx")
path_LIFEHeart_GE = paste0(basicpath,"/02_qcdaten/1906_LifeB3_gx/ExpressionSet_Prepro_v3_LIFE_Heart_filtered_IPA.RData")

# genetic data
path_LIFEAdult_Genetics = paste0(basicpath,"/02_projekte/1612_lifea1_genotypisierungsrunde_3/Imputation/08_pgen_TOPMed/results/LIFE-Adult_TOPMed_minimal_filtered_chr1to23")
path_LIFEHeart_Genetics = paste0(basicpath,"/02_projekte/1509_lifeb3_axiom/10_pgen_TOPMed/results/LIFE-Heart_TOPMed_minimal_filtered_chr1to23")
path_LIFEAdult_GeneticPCs = paste0(basicpath,"/02_qcdaten/LifeA1_Rd3/SamplesAsImputed.RData")
path_LIFEHeart_GeneticPCs = paste0(basicpath,"/02_qcdaten/LifeB3_Rd2/SamplesAsImputedB3.RData")

# results from data preparation
path_LIFEprepped = paste0(projectpath,"_data/dataQC/")

#############################
# path to summary statistics 
#############################
# Steroid hormone GWAS Pott et al. 
path_SumStats_2021 = paste0(basicpath,"/02_projekte/1910_lifea1b3_sasha/2004_Meta_GWAS/SumStats/")
path_SumStats_2019 = paste0(basicpath,"/02_projekte/1604_lifeb3_gwas_sexhormone/1807_meta_gwas/SummaryStats/")

path_Ruth_TT_2020 = paste0(projectpath,"_data/SummaryStatistics/32042192")
path_Schmitz_E2_2021 = paste0(projectpath,"_data/SummaryStatistics/34255042")
path_Chan_CORT_2024 = paste0(projectpath,"_data/SummaryStatistics/38525495")

path_SumStats_QC = paste0(projectpath,"_data/SummaryStatistics_QC/")

#############################
# path to other data sets 
#############################
path_PathwayGenesKEGG = paste0(projectpath,"_data/otherData/230702_KEGG_genes_hsa00140.txt")
