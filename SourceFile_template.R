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

suppressPackageStartupMessages(library(ivreg))
suppressPackageStartupMessages(library(pgenlibr))
suppressPackageStartupMessages(library(MendelianRandomization))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

basicpath = "/Path/To/Projects/"
projectpath = "/Path/To/This/Proejct/"

#############################
# Other Tools 
#############################
# PLINK2 path 
path_plink2 = "/Path/To/plink2.0/20210203/unix_64/plink2" 
path_plink1.9 = "/Path/To/plink1.9/vs210606_stable_beta_6_24/unix64/plink"

#############################
# Other helper functions 
#############################
# please check with Holger Kirsten for access (toolbox via github, limma function per request)
suppressPackageStartupMessages(library(toolboxH))
func_TWAS = paste0(basicpath,"LimmaFunktion_HolgerKirsten.R") 

# these are my functions
func_TWASExtract = "helperfunctions/myExtractionFunction.R"

#############################
# path to LIFE data 
#############################
# all individual level data is not part of the repository and cannot be shared!
# Please contact LIFE data management if you want access to this data. 
# Data access via the LIFE project agreement 785
#
# phenotype data
path_LIFEAdult = paste0(projectpath,"_data/pv785/data_adult/")
path_LIFEHeart = paste0(projectpath,"_data/pv785/data_heart/")
path_LIFEHeart_old = paste0(projectpath,"_data/pv505/data_heart/")

# gene expression data
path_LIFEAdult_GE = paste0(basicpath,"/02_qcdaten/LIFE_Adult/ExpressionSet.RData")
path_LIFEAdult_GE_doku = paste0(basicpath,"/02_qcdaten/LIFE_Adult/DokuTab.xlsx")
path_LIFEHeart_GE = paste0(basicpath,"/02_qcdaten/LIFE_HEART/ExpressionSet.RData")

# genetic data
path_LIFEAdult_Genetics = paste0(basicpath,"/02_qcdaten/LIFE_Adult/TOPMed")
path_LIFEHeart_Genetics = paste0(basicpath,"/02_qcdaten/LIFE_Heart/TOPMed")
path_LIFEAdult_GeneticPCs = paste0(basicpath,"/02_qcdaten/LIFE_Adult/SamplesAsImputed.RData")
path_LIFEHeart_GeneticPCs = paste0(basicpath,"/02_qcdaten/LIFE_Heart/SamplesAsImputedB3.RData")

# results from data preparation
path_LIFEprepped = paste0(projectpath,"_data/dataQC/")

#############################
# path to UKB data 
#############################
# all individual level data is not part of the repository and cannot be shared!
# Data accesss via UKB application number 98032 
# 
UKB_SNP_data = "PATH/TO/UKB/genotypes-imputed/"
UKB_phenotypes = "PATH/TO/UKBB/phenotypes/ukb672224.tab"
UKB_phenotypes_filtered = "PATH/TO/UKB_projectFiltered_phenotypes/"
UKB_genotypes_filtered = "PATH/TO/UKB_projectFiltered_genotypes/"
