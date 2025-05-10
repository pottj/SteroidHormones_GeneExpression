#' ---
#' title: "PGS - Calculating PGS - CORT special"
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
#' CORT is monogenetic, and trying to use a genome-wide or pathway-wide score resulted in no significant results. Hence, I want to use the Cortisol Binding Globulin (CBG, aka *SERPINA6*) region only, and test if this helps. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Restrict to CBG ####
#' ***
files = list.files(path = path_SumStats_QC)
files = files[grepl(".gz",files)]
files = files[grepl("Pathway",files)]
files = files[grepl("CORT",files)]

CORT = fread(paste0(path_SumStats_QC,files))
CORT = CORT[genes=="SERPINA6"]

dummy = gsub(".txt.gz","",files)
dummy = gsub("PathwayRegions","SERPINA6",dummy)
outfile = paste0(path_SumStats_QC,dummy,".txt")
fwrite(CORT, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F, na = NA, dec = ".")
gzip(outfile,destname = paste0(outfile, ".gz"))

write.table(CORT$ID,file = "../results/02_PGS_07_SNPList_raw.txt", 
            col.names = F, row.names = F, quote = F)

#' # Create input parameters ####
#' ***
#' Get summary statistics file names and study IDs
studies = c("Adult","Heart")

#' Create a file containing the different p-value thresholds for inclusion of SNPs in the PGS.
thresholds = c(1e-3,1e-4,1e-5,1e-6)

#' # Loop ####
#' ***
for(i in 1:length(files)){
  #i=1
  dummy = gsub(".txt.gz","",files[i])
  dummy = gsub("PathwayRegions","SERPINA6",dummy)
  hormone = gsub("_.*","",dummy)
  setting = gsub("_SERPINA6.*","",dummy)
  setting = gsub(".*_","",setting)
  message("Working on ",hormone, " (setting: ",setting,")")
  
  # Step 0: extract the SNP ID (column 14) and p-values (column 7) from the summary statistics file
  myCall0 = paste0("zcat ",path_SumStats_QC,dummy,".txt.gz | awk '{print $14,$7}' > ",path_SumStats_QC,hormone,"_",setting,"_SERPINA6.pvalue")
  system(myCall0)
  
  for(j in 1:length(studies)){
    #j=1
    time1<-Sys.time()
    study = studies[j]
    message("Working on ",hormone, " (setting: ",setting,") in the LIFE-",study," study")
    
    if(study == "Adult"){
      path_genetics = path_LIFEAdult_Genetics
    }else{
      path_genetics = path_LIFEHeart_Genetics
    }
    
    # Step 1: make bed file for clumping
    myCall1 = paste0(path_plink2,
                     " --pfile ",path_genetics,
                     " --extract ../results/02_PGS_07_SNPList_raw.txt",
                     " --maf 0.01",
                     " --hwe 1e-6",
                     " --geno 0.01",
                     " --mind 0.01", 
                     " --mach-r2-filter 0.8 2",
                     " --make-bed",
                     " --out ",path_LIFEprepped,"genetics/",study,"_QC_SERPINA6")
    myCall1
    system(myCall1)
    
    # Step 2: Clumping (PLINK 1.9)
    if(study == "Adult"){
      myCall2 = paste0(path_plink1.9,
                       " --bfile ",path_LIFEprepped,"genetics/",study,"_QC_SERPINA6",
                       " --clump-p1 1",
                       " --clump-r2 0.1",
                       " --clump-kb 250",
                       " --clump ",path_SumStats_QC,files[i],
                       " --clump-snp-field ID",
                       " --clump-field P",
                       " --out ../results/02_PGS_07_PLINK_PGS_SERPINA6/",study,"_Clumping_",hormone,"_",setting)
      myCall2
      system(myCall2)
    }
    
    # Step 3: extract clumped SNPs
    clumpedSNPs = fread(paste0("../results/02_PGS_07_PLINK_PGS_SERPINA6/Adult_Clumping_",hormone,"_",setting,".clumped"),header = T)
    CORT2 = copy(CORT)
    CORT2 = CORT2[ID %in% clumpedSNPs$SNP]
    setorder(CORT2,BP)
    
    # Step 4: load pgen file into R
    pvar1 = NewPvar(paste0(path_genetics, '.pvar'))
    pgen = NewPgen(paste0(path_genetics,'.pgen'), pvar=pvar1)
    pvar = fread(paste0(path_genetics,'.pvar'))
    psam = fread(paste0(path_genetics,'.psam'))
    setnames(pvar,"#CHROM","CHR")
    
    filt = is.element(pvar$ID,CORT2[P<=0.001, ID])
    table(filt)
    myNRs = 1:dim(pvar)[1]
    myNRs = myNRs[filt]
    geno_mat <- ReadList(pgen, variant_subset = myNRs , meanimpute=F)
    dim(geno_mat)
    colnames(geno_mat) = pvar[filt,ID]
    rownames(geno_mat) = psam$IID
    
    # Step 5: check alleles real quick again
    pvar = pvar[filt,]
    CORT2 = CORT2[P<0.001,]
    table(pvar$ID %in% CORT2$ID)
    stopifnot(pvar$ID == CORT2$ID)
    stopifnot(pvar$ALT == CORT2$EA)
    
    # Step 6: calculate scores
    score_mat = matrix(nrow = dim(geno_mat)[1],ncol = length(thresholds))
    for(k in 1:length(thresholds)){
      #k=1
      filt = CORT2$P<thresholds[k]
      table(filt)
      weights = CORT2$BETA[filt]
      score = geno_mat[,filt] %*% weights
      score_mat[,k] = score[,1]
    }
    
    # Step 7: merge to samples as imputed
    if(studies[j]=="Adult"){
      loaded = load(path_LIFEAdult_GeneticPCs)
      SamplesAsImputed = get(loaded)
      SamplesAsImputed = SamplesAsImputed[filterGonosom==F,]
      SamplesAsImputed = SamplesAsImputed[,c(2,5,8:17)]
      #names(SamplesAsImputed)
    }else{
      loaded = load(path_LIFEHeart_GeneticPCs)
      SamplesAsImputed = get(loaded)
      setDT(SamplesAsImputed)
      SamplesAsImputed = SamplesAsImputed[filterGonosom==F,]
      SamplesAsImputed = SamplesAsImputed[,c(6,10,12:21)]
      setnames(SamplesAsImputed,"Aliquot_Genetik","Aliquot")
      #names(SamplesAsImputed)
    }
    
    matched = match(SamplesAsImputed$Aliquot,psam$IID) 
    table(SamplesAsImputed$Aliquot == psam$IID[matched])
    SamplesAsImputed = cbind(SamplesAsImputed,score_mat[matched,])
    names(SamplesAsImputed)[13:16] = paste0("CORT_combined_1e-",c(3:6))
    
    # Step 8: save object 
    save(SamplesAsImputed, file = paste0(path_LIFEprepped,"genetics/",studies[j],"_PGS_SERPINA6.RData"))
    
    # Step 9: check time
    time_dif = round(difftime(Sys.time(),time1,units = "mins"),3)
    print(time_dif)
  }
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
