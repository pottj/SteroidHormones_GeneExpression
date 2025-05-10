#' ---
#' title: "Get UKB data - SNP selection"
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
library(foreach)

#' # Parameter settings ####
#' ***
UKB_SNP_data = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/"
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"

#' # Save SNPs #### 
#' ***
#' Please check project SteroidHormones_GeneExpression for the selection process! 
#' 
#' The selected SNPs are 
#' 
#' - genome-wide significantly associated with testosterone levels in men in Ruth et al. 
#' - within/near genes of the steroid hormone biosynthesis pathway
#' 
load("../temp/proteomics_instruments.RData")
write.table(unique(erg$rsID),
            file = paste0(data_QC,"/01_Prep_02_SNPList.txt"), 
            col.names = F, row.names = F, quote = F)

#' # Create PLINK2 calls ####
#' ***
#' I generate the calls here, and then copy-paste them into a slurm script. 
#' 
#' 1) create temp folder
#' 2) create bgen per chromosome (stored in temp)
#' 3) merge bgen files into one file (stored in temp)
#' 4) create pgen merged file (stored **not** in temp)
#' 5) remove temp folder (do not store all the chromosome-wise data)
#' 
call1 = paste0("mkdir ",data_QC,"/temp")
print(call1)

myCHR = unique(erg$CHR)
myCHR = myCHR[!is.na(myCHR)]
myCHR = myCHR[order(myCHR)]

dumTab = foreach(i = 1:length(myCHR))%do%{
  #i=1
  call2 = paste0("plink2", 
                 " --bgen ",UKB_SNP_data,"/ukb22828_c",myCHR[i],"_b0_v3.bgen",
                 " 'ref-last'", 
                 " --sample ",UKB_SNP_data,"ukb22828_c",myCHR[i],"_b0_v3_s487160.sample",
                 " --chr ", myCHR[i],
                 " --keep-fam ",data_QC,"/01_Prep_01_SampleList.txt",
                 " --extract ",data_QC,"/01_Prep_02_SNPList.txt", 
                 " --mach-r2-filter 0.8 2",
                 " --maf 0.01", 
                 " --threads 20",
                 " --export bgen-1.2 bits=8 id-delim='-'",
                 " --out ",data_QC,"/temp/UKB_Testo_chr",myCHR[i])
  print(call2)
  out_filename = paste0(data_QC,"/temp/UKB_Testo_chr",myCHR[i],".bgen")
  out_filename
}

call3 = c("cat-bgen -clobber -g")

for(i in 1:length(myCHR)){
  #i=1
  call3 = paste(call3,dumTab[[i]])
}

call3 = paste0(call3, " -og ",data_QC, "/temp/UKB_Testo_merged.bgen")
print(call3)

call4 = paste0("plink2",
               " --bgen ",data_QC, "/temp/UKB_Testo_merged.bgen",
               " 'ref-last'",
               " --sample ",data_QC, "/temp/UKB_Testo_chr2.sample",
               " --make-pgen --out ",data_QC,"/01_Prep_02_UKB_Testo_merged")
print(call4)

call5 = paste0("rm -rf ",data_QC,"/temp")
print(call5)

call6 = paste0("plink2",
               " --pfile ",data_QC, "/01_Prep_02_UKB_Testo_merged",
               " --freq --out ",data_QC,"/01_Prep_02_UKB_Testo_merged_AF")
print(call6)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
