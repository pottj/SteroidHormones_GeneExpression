#' ---
#' title: "PGS - Calculating PGS"
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
#' Here, I need to create a loop to do for each hormone and study
#' 
#' - clumping (PLINK 1.9)
#' - score calculation (PLINK 2.0)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Create input parameters ####
#' ***
#' Get summary statistics file names and study IDs
files = list.files(path = path_SumStats_QC)
files = files[grepl(".gz",files)]
studies = c("Adult","Heart")

myTab = data.table(file = rep(files,each=2),
                   study = rep(studies,5))
myTab[,hormone := gsub("_.*","",file)]
myTab[,setting := gsub("[.].*","",file)]
myTab[,setting := gsub(".*_","",setting)]

#' Create a file containing the different p-value thresholds for inclusion of SNPs in the PGS.
dummy = c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
if(file.exists("../results/03_range_list")){
  # I might end up rerunning this with different thresholds
  myCall = paste0("rm ../results/03_range_list")
  system(myCall)
}
for(i in 1:length(dummy)){
  myCall = paste0("echo '",dummy[i],' 0 ',dummy[i],"' >> ../results/03_range_list")
  system(myCall)
}

#' # Loop ####
#' ***
for(i in 1:length(files)){
  #i=1
  hormone = gsub("_.*","",files[i])
  setting = gsub("[.].*","",files[i])
  setting = gsub(".*_","",setting)
  message("Working on ",hormone, " (setting: ",setting,")")
  
  # Step 0: extract the SNP ID (column 14) and p-values (column 7) from the summary statistics file
  myCall0 = paste0("zcat ",path_SumStats_QC,files[i]," | awk '{print $14,$7}' > ",path_SumStats_QC,hormone,"_",setting,".pvalue")
  system(myCall0)
  
  for(j in 1:length(studies)){
    #j=1
    time1<-Sys.time()
    study = studies[j]
    message("Working on ",hormone, " (setting: ",setting,") in the LIFE-",study," study")
    
    # Step 1: Clumping (PLINK 1.9)
    myCall1 = paste0(path_plink1.9,
                     " --bfile ",path_LIFEprepped,"genetics/",study,"_QC",
                     " --clump-p1 1",
                     " --clump-r2 0.1",
                     " --clump-kb 250",
                     " --clump ",path_SumStats_QC,files[i],
                     " --clump-snp-field ID",
                     " --clump-field P",
                     " --out ../results/03_PLINK_PGS/",study,"_Clumping_",hormone,"_",setting)
    myCall1
    system(myCall1)
    
    # Step 2: extract SNPs out of result file (these are the clumped independent SNPs)
    myCall2 = paste0("awk 'NR!=1{print $3}' ../results/03_PLINK_PGS/",study,"_Clumping_",hormone,"_",setting,".clumped >  ../results/03_PLINK_PGS/",study,"_Clumping_",hormone,"_",setting,".valid.snp")
    system(myCall2)
   
    # Step 3: get scores using the ID (col 14), effect allele (col 4) and effect size (col 8) of the summary statistics file
    if(study == "Adult"){
      path_genetics = path_LIFEAdult_Genetics
    }else{
      path_genetics = path_LIFEHeart_Genetics
    }
    myCall3 = paste0(path_plink2,
                     " --pfile ",path_genetics,
                     " --extract ../results/03_PLINK_PGS/",study,"_Clumping_",hormone,"_",setting,".valid.snp",
                     " --score ",path_SumStats_QC,files[i]," 14 4 8 header",
                     " --q-score-range ../results/03_range_list ",path_SumStats_QC,hormone,"_",setting,".pvalue",
                     " --out ../results/03_PLINK_PGS/",study,"_Score_",hormone,"_",setting)
    myCall3
    system(myCall3)
    
    # Step 4: check time
    time_dif = round(difftime(Sys.time(),time1,units = "mins"),3)
    myTab[file == files[i] & study == studies[j], time := time_dif]
  }
}
myTab

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
