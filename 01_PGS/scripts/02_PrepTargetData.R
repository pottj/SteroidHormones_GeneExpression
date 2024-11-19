#' ---
#' title: "PGS - QC of target data"
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
#' Here, I check the genetic data of LIFE-Adult and LIFE-Heart before getting the genome-wide PGS per hormone: 
#' 
#' - Standard GWAS QC: MAF>1%, hwe, imputation info, missingness
#' - LD pruning
#' - save as bed file (necessary for clumping)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Adult ####
#' ***
#' Step 1: standard QC

myCall1 = paste0(path_plink2,
                 " --pfile ",path_LIFEAdult_Genetics,
                 " --extract ../results/01_SNPList_raw.txt",
                 " --maf 0.01",
                 " --hwe 1e-6",
                 " --geno 0.01",
                 " --mind 0.01", 
                 " --mach-r2-filter 0.8 2",
                 " --write-snplist", 
                 " --make-just-fam", 
                 " --out ../results/02_PLINK_prep/Adult_QC_SNPList")
myCall1
system(myCall1)

#' Step 2: pruning

myCall2 = paste0(path_plink2,
                 " --pfile ",path_LIFEAdult_Genetics,
                 " --extract ../results/02_PLINK_prep/Adult_QC_SNPList.snplist",
                 " --keep ../results/02_PLINK_prep/Adult_QC_SNPList.fam", 
                 " --indep-pairwise 200 50 0.25", 
                 " --out ../results/02_PLINK_prep/Adult_QC_Pruning")
myCall2
system(myCall2)

#' Step 3: make bed-file for later use (stored in dataQC, not tracked)
#' 
myCall3 = paste0(path_plink2,
                 " --pfile ",path_LIFEAdult_Genetics,
                 " --extract ../results/02_PLINK_prep/Adult_QC_Pruning.prune.in",
                 " --make-bed",
                 " --out ",path_LIFEprepped,"genetics/Adult_QC")
myCall3
system(myCall3)


#' # Heart ####
#' ***
#' Note: For Heart, I only generate the bad files, as I want to use the same SNPs for the score as in Adult. Study-specific pruning might prevent that from happening!
#' 
#' Step 1: standard QC

# myCall4 = paste0(path_plink2,
#                  " --pfile ",path_LIFEHeart_Genetics,
#                  " --extract ../results/01_SNPList_raw.txt",
#                  " --maf 0.01",
#                  " --hwe 1e-6",
#                  " --geno 0.01",
#                  " --mind 0.01", 
#                  " --mach-r2-filter 0.8 2",
#                  " --write-snplist", 
#                  " --make-just-fam", 
#                  " --out ../results/02_PLINK_prep/Heart_QC_SNPList")
# myCall4
# system(myCall4)

#' Step 2: pruning

# myCall5 = paste0(path_plink2,
#                  " --pfile ",path_LIFEHeart_Genetics,
#                  " --extract ../results/02_PLINK_prep/Heart_QC_SNPList.snplist",
#                  " --keep ../results/02_PLINK_prep/Heart_QC_SNPList.fam", 
#                  " --indep-pairwise 200 50 0.25", 
#                  " --out ../results/02_PLINK_prep/Heart_QC_Pruning")
# myCall5
# system(myCall5)

#' Step 3: make bed-file for later use (stored in dataQC, not tracked)
#' 
myCall6 = paste0(path_plink2,
                 " --pfile ",path_LIFEHeart_Genetics,
                 # " --extract ../results/02_PLINK_prep/Heart_QC_Pruning.prune.in",
                 " --extract ../results/02_PLINK_prep/Adult_QC_Pruning.prune.in",
                 " --make-bed",
                 " --out ",path_LIFEprepped,"genetics/Heart_QC")
myCall6
system(myCall6)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
