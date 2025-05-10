#' ---
#' title: "PGS - Calculating PGS restricted to pathway region"
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
#' First, I have to load the summary statistic files and restrict them to the relevant gene regions. Then I repeat the steps of script 03 to calculate the PGS. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Restrict to pathway regions ####
#' ***
files = list.files(path = path_SumStats_QC)
files = files[grepl(".gz",files)]
files = files[!grepl("PathwayRegions",files)]
files = files[!grepl("SERPINA6",files)]

load("../../00_dataPrep/results/01_prep_03_PathwayGenes_byCytoband.RData")

for(i in 1:length(files)){
  #i=1
  message("Working on file ",files[i]," (number ",i," of ",length(files),")")
  
  # Step 1: load summary statistics
  erg1 = fread(paste0(path_SumStats_QC,files[i]))
  
  # Step 2: all SNPs of the pathway gene regions
  dumTab2 = foreach(j = 1:dim(PathwayGenes_filt)[1])%do%{
    #j=1
    message("       Working on locus number ",j)
    erg4 = copy(erg1)
    myLocus = copy(PathwayGenes_filt)
    myLocus = myLocus[j,]
    
    erg4 = erg4[CHR == myLocus$CHR,]
    erg4 = erg4[BP >= myLocus$region_start_collapsed,]
    erg4 = erg4[BP <= myLocus$region_end_collapsed,]
    
    erg4[,cyto:= myLocus$cytoband]
    erg4[,genes:= myLocus$genes]
    
    erg4
  }
  erg3 = rbindlist(dumTab2)
  
  # Step 3: save as txt file and zip it
  dummy = gsub(".txt.gz","",files[i])
  outfile = paste0(path_SumStats_QC,dummy,"_PathwayRegions.txt")
  fwrite(erg3, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F, na = NA, dec = ".")
  gzip(outfile,destname = paste0(outfile, ".gz"))
  
}

#' # Create input parameters ####
#' ***
#' Get summary statistics file names and study IDs
files = list.files(path = path_SumStats_QC)
files = files[grepl(".gz",files)]
files = files[grepl("Pathway",files)]
studies = c("Adult","Heart")

myTab = data.table(file = rep(files,each=2),
                   study = rep(studies,5))
myTab[,hormone := gsub("_.*","",file)]
myTab[,setting := gsub("[.].*","",file)]
myTab[,setting := gsub(".*_","",setting)]

#' Create a file containing the different p-value thresholds for inclusion of SNPs in the PGS.
dummy = c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
if(file.exists("../results/02_PGS_05_range_list")){
  # I might end up rerunning this with different thresholds
  myCall = paste0("rm ../results/02_PGS_05_range_list")
  system(myCall)
}
for(i in 1:length(dummy)){
  myCall = paste0("echo '",dummy[i],' 0 ',dummy[i],"' >> ../results/02_PGS_05_range_list")
  system(myCall)
}

#' # Loop ####
#' ***
for(i in 1:length(files)){
  #i=1
  hormone = gsub("_.*","",files[i])
  setting = gsub("_Pathway.*","",files[i])
  setting = gsub(".*_","",setting)
  message("Working on ",hormone, " (setting: ",setting,")")
  
  # Step 0: extract the SNP ID (column 14) and p-values (column 7) from the summary statistics file
  myCall0 = paste0("zcat ",path_SumStats_QC,files[i]," | awk '{print $14,$7}' > ",path_SumStats_QC,hormone,"_",setting,"_Pathway.pvalue")
  system(myCall0)
  
  # Step 1: Clumping (PLINK 1.9) using LIFE-Adult data
  myCall1 = paste0(path_plink1.9,
                   " --bfile ",path_LIFEprepped,"genetics/Adult_QC",
                   " --clump-p1 1",
                   " --clump-r2 0.1",
                   " --clump-kb 250",
                   " --clump ",path_SumStats_QC,files[i],
                   " --clump-snp-field ID",
                   " --clump-field P",
                   " --out ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Clumping_",hormone,"_",setting)
  myCall1
  system(myCall1)
  
  # Step 2: extract SNPs out of result file (these are the clumped independent SNPs)
  myCall2 = paste0("awk 'NR!=1{print $3}' ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Clumping_",hormone,"_",setting,".clumped >  ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Clumping_",hormone,"_",setting,".valid.snp")
  system(myCall2)
  
  # Step 3: get scores using the ID (col 14), effect allele (col 4) and effect size (col 8) of the summary statistics file
  
  myCall3a = paste0(path_plink2,
                    " --pfile ",path_LIFEAdult_Genetics,
                    " --extract ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Clumping_",hormone,"_",setting,".valid.snp",
                    " --score ",path_SumStats_QC,files[i]," 14 4 8 header",
                    " --q-score-range ../results/02_PGS_03_range_list ",path_SumStats_QC,hormone,"_",setting,".pvalue",
                    " --out ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Score_",hormone,"_",setting)
  myCall3a
  system(myCall3a)
  
  myCall3b = paste0(path_plink2,
                    " --pfile ",path_LIFEHeart_Genetics,
                    " --extract ../results/02_PGS_05_PLINK_PGS_pathway/Adult_Clumping_",hormone,"_",setting,".valid.snp",
                    " --score ",path_SumStats_QC,files[i]," 14 4 8 header",
                    " --q-score-range ../results/02_PGS_03_range_list ",path_SumStats_QC,hormone,"_",setting,".pvalue",
                    " --out ../results/02_PGS_05_PLINK_PGS_pathway/Heart_Score_",hormone,"_",setting)
  myCall3b
  system(myCall3b)
  
}
myTab

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
