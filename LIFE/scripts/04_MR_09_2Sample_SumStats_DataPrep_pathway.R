#' ---
#' title: "MR - 2-sample - summary statistics"
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
#' Plan: 
#' 
#' - get instruments from summary statistics
#' - load eqtlgen data set
#' - check for overlapping SNP IDs
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

source("../../helperfunctions/getSmallestDist.R")

#' # Load eQTLGen data ####
#' ***
eqtlgen = fread("../../../_data/SummaryStatistics/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz")

all_genes = eqtlgen[,.N,by=c("Gene","GeneSymbol")]
all_SNPs = eqtlgen[,.N,by = c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele")]

save(all_SNPs,file = "../temp/04_MR_eQTLGen_allSNPs.RData")
save(all_genes,file = "../temp/04_MR_eQTLGen_allgenes.RData")

load("../temp/04_MR_eQTLGen_allSNPs.RData")
load("../temp/04_MR_eQTLGen_allgenes.RData")
load("../results/01_prep_03_PathwayGenes_byCytoband.RData")

#' # Get instruments ####
#' ***
myFiles = list.files(path = path_SumStats_QC)
myFiles = myFiles[!grepl("athway",myFiles)]
myFiles = myFiles[!grepl("SERPINA",myFiles)]
myFiles = myFiles[!grepl("pval",myFiles)]
myFiles = myFiles[!grepl("genomewide",myFiles)]

dumTab1 = foreach(i = 1:5)%do%{
  #i=1
  erg2 = fread(paste0(path_SumStats_QC,myFiles[i]))
  
  erg2[,Zscore := BETA/SE]
  erg2[,Zscore_abs := abs(BETA/SE)]
  table(is.element(erg2$rsID,all_SNPs$SNP),erg2$Zscore_abs>abs(qnorm(1e-6)))
  
  erg3 = erg2[rsID %in% all_SNPs$SNP,]
  #erg3 = erg3[Zscore_abs>abs(qnorm(1e-6)),]
  
  # Step 2: all SNPs of the pathway gene regions
  dumTab2 = foreach(j = 1:dim(PathwayGenes_filt)[1])%do%{
    #j=1
    message("       Working on locus number ",j)
    erg4 = copy(erg3)
    myLocus = copy(PathwayGenes_filt)
    myLocus = myLocus[j,]
    
    erg4 = erg4[CHR == myLocus$CHR,]
    erg4 = erg4[BP >= myLocus$region_start_collapsed,]
    erg4 = erg4[BP <= myLocus$region_end_collapsed,]
    
    erg4[,cyto:= myLocus$cytoband]
    erg4[,genes:= myLocus$genes]
    
    erg4
  }
  erg5 = rbindlist(dumTab2)
  
  if(dim(erg5)[1]>0){
    # priority pruning
    setorder(erg5,CHR,BP)
    myCHRs = unique(erg5$CHR)
    
    result.22 = foreach(s2 = myCHRs) %do% {
      # s2 = myCHRs[1]
      subdata2 = copy(erg5)
      subdata2 = subdata2[CHR == s2, ]
      setkey(subdata2, BP)
      
      if(dim(subdata2)[1]<=1){
        subdata2[, keep := T]
        subdata2[, NR_SNPs := 0]
      }else{
        subdata2[, keep := NA]
        subdata2[, NR_SNPs := as.numeric(NA)]
        
        smallestDist = getSmallestDist(subdata2[, BP])
        while(smallestDist < 1000000) {
          #minP = min(subdata2[is.na(keep), p])
          #maxLogP = max(subdata2[is.na(keep), pvalue_neg_log10_GC])
          maxAbsZscore = max(subdata2[is.na(keep), Zscore_abs])
          myPOS = subdata2[maxAbsZscore == Zscore_abs & is.na(keep), BP]
          if(length(myPOS)>1){
            myPOS = myPOS[1]
          }
          subdata2[BP == myPOS, keep := T]
          
          #filter for SNPs that can stay within the set (outside the +- 1MB range or keep==T)
          myFilt = (subdata2[, BP] < (myPOS - 1000000)) | 
            (subdata2[, BP] > (myPOS + 1000000)) | 
            subdata2[, keep] 
          myFilt[is.na(myFilt)] = FALSE
          subdata2 = subdata2[myFilt == TRUE, ]
          
          subdata2[BP == myPOS, NR_SNPs := sum(myFilt==F)]
          smallestDist = getSmallestDist(subdata2[, BP])
        }
        
        #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
        subdata2[is.na(keep), NR_SNPs := 0]
        subdata2[is.na(keep), keep := TRUE]
      }
      
      subdata2
    }
    erg6 = rbindlist(result.22)
    myFile2 = gsub(".txt.gz","",myFiles[i])
    outfile = paste0(path_SumStats_QC,myFile2,"_pathway_eqtlgenMatched.txt")
    fwrite(erg6, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F, na = NA, dec = ".")
  }
  erg6
}
erg = rbindlist(dumTab1,use.names = T,fill=T)
erg[,table(Zscore_abs>5,phenotype,setting)]

#' Summary: MR can only be done for TT in both sexes, and E2 in men
erg[,trait := paste(phenotype,setting,sep="_")]
erg = erg[trait %in% c("E2_men","TT_men","TT_women")]
erg = erg[Zscore_abs>5,]
save(erg,file = "../temp/04_MR_UKB_IVs_pathway.RData")

table(is.element(erg$rsID,all_SNPs$SNP))

#' # Prep eQTL data ####
#' ***
#' Restrict to instruments, check alleles, get LIFE-Adult EAF, calculate beta and SE 
#' 
filtered_SNPs = copy(all_SNPs)
filtered_SNPs = filtered_SNPs[SNP %in% erg$rsID,]

load("../results/02_PGS_01_LIFE_pvar_filtered.RData")
table(is.element(erg$ID,pvar$ID))
pvar = pvar[ID %in% erg$ID,]
matched1 = match(pvar$ID,erg$ID)
stopifnot(pvar$ID == erg$ID[matched1])
pvar[,rsID := erg[matched1,rsID]]
stopifnot(pvar$rsID %in% filtered_SNPs$SNP)

matched2 = match(pvar$rsID,filtered_SNPs$SNP)
filtered_SNPs = filtered_SNPs[matched2,]

stopifnot(filtered_SNPs$SNP == pvar$rsID)
stopifnot(filtered_SNPs$SNPChr == pvar$CHR)
plot(filtered_SNPs$SNPPos, pvar$POS)

table(filtered_SNPs$AssessedAllele == pvar$ALT,filtered_SNPs$OtherAllele == pvar$REF)
table(filtered_SNPs$OtherAllele == pvar$ALT,filtered_SNPs$AssessedAllele == pvar$REF)
table(erg[matched1,EA] == pvar$ALT,erg[matched1,OA] == pvar$REF)
filtered_SNPs[,switchAlleles := F]
filtered_SNPs[OtherAllele == pvar$ALT, switchAlleles := T]
filtered_SNPs[,EAF_LIFE := pvar$AF_Adult]
filtered_SNPs[,MAF_LIFE := EAF_LIFE]
filtered_SNPs[EAF_LIFE>0.5,MAF_LIFE := 1-EAF_LIFE]
plot(filtered_SNPs$EAF_LIFE,filtered_SNPs$MAF_LIFE)

eqtlgen_filtered = eqtlgen[SNP %in% filtered_SNPs$SNP,]

dumTab2 = foreach(i = 1:dim(filtered_SNPs)[1])%do%{
  #i=1
  mySNP = filtered_SNPs[i,]
  erg5 = copy(eqtlgen_filtered)
  erg5 = erg5[SNP == mySNP$SNP,]
  
  if(mySNP$switchAlleles==T){
    erg5[,EA := OtherAllele]
    erg5[,OA := AssessedAllele]
    erg5[,Zscore_harmonized := Zscore * (-1)]
  }else{
    erg5[,EA := AssessedAllele]
    erg5[,OA := OtherAllele]
    erg5[,Zscore_harmonized := Zscore]
  }
  
  erg5[, beta := Zscore_harmonized / sqrt(NrSamples * 2 * mySNP$MAF_LIFE * (1-mySNP$MAF_LIFE))]
  erg5[, SE := beta/Zscore_harmonized]
  erg5 = erg5[beta!=0,]
  erg5
}
eqtlgen_harmonized = rbindlist(dumTab2)
save(eqtlgen_harmonized,file = "../temp/04_MR_eqtlgen_harmonized_pathway.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

