#' ---
#' title: "Get CORT SNPs"
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
#' In this script, I will load the cortisol data from Chan and Wu (2024), and filter for 
#' 
#' - the gene regions defined in a previous script, and
#' - availability in LIFE (TopMed imputation).
#' 
#' Then, I will extract per sex-setting the best-associated SNP per cytoband. This gives me a list of 2 x 37 SNPs (not necessarily the same per setting) 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_angmar.R")
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP150.GRCh38))
.libPaths()

#' # Load data ####
#' ***
load("../results/01_PathwayGenes_byCytoband.RData")

patterns = gsub(".*/","",path_Chan_CORT_2024)
path_Chan_CORT_2024 = gsub(patterns,"",path_Chan_CORT_2024)
statistics = list.files(path = path_Chan_CORT_2024, pattern = patterns)

ToDoList = data.table(filenames = statistics,
                      hormone = "CORT",
                      setting = c("sex-combined"),
                      sampleSize = c(32981))

myNames_old = c("hormone","setting",
                "rsid","chr","pos","A1","A2","Freq1",
                "beta","se","P","zscore","N")
myNames_new = c("hormone","setting",
                "rsID","CHR","POS","EA","OA","EAF",
                "beta","SE","pvalue","zscore","sampleSize")

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myTime = Sys.time()
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$hormone," in ",myRow$setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(path_Chan_CORT_2024,myRow$filenames)
  erg1 = fread(myfn1)
  
  erg1[,hormone := myRow$hormone]
  erg1[,setting := myRow$setting]
  erg1[,zscore := beta/se]
  
  stopifnot(myNames_old %in% names(erg1))
  colsOut<-setdiff(colnames(erg1),myNames_old)
  erg1[,get("colsOut"):=NULL]
  setcolorder(erg1,myNames_old)
  names(erg1) = myNames_new
  
  # lift over from hg19 to hg38!
  erg2 = copy(erg1)
  erg2 = erg2[grepl("rs",rsID),]
  mySNPs = erg2$rsID
  snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
  myLift = snpsById(snps, mySNPs, ifnotfound="drop")
  
  myPos = pos(myLift)
  myIDs = myLift$RefSNP_id
  
  matched = match(mySNPs,myIDs)
  erg2[,pos_b38 := myPos[matched]]
  erg2 = erg2[!is.na(pos_b38)]
  erg2[,POS := NULL]
  setnames(erg2,"pos_b38","POS")
  setcolorder(erg1,myNames_new)
  
  # get list 3: all SNPs of the pathway gene regions
  dumTab2 = foreach(j = 1:dim(PathwayGenes_filt)[1])%do%{
    #j=1
    message("Working on locus number ",j)
    erg4 = copy(erg2)
    myLocus = copy(PathwayGenes_filt)
    myLocus = myLocus[j,]
    
    erg4 = erg4[CHR == myLocus$CHR,]
    erg4 = erg4[POS >= myLocus$region_start_collapsed,]
    erg4 = erg4[POS <= myLocus$region_end_collapsed,]
    
    erg4[,cyto:= myLocus$cytoband]
    erg4[,genes:= myLocus$genes]
    
    erg4
  }
  erg3 = rbindlist(dumTab2)
  erg3
  
}
myTab_complete = rbindlist(dumTab1)
myTab_complete[,MAF := EAF]
myTab_complete[EAF>0.5,MAF := 1-EAF]

save(myTab_complete,file ="../temp/04_CORT_unfiltered.RData")

#' # Check LIFE-Adult ####
#' ***
#' I will load the pvar file and check the overlap. 
#' 
pvar = fread(paste0(path_LIFEAdult_Genetics,".pvar"))
head(pvar)
myTab_complete[,dumID1 := paste0("chr",CHR,":",POS,":",EA,":",OA)]
myTab_complete[,dumID2 := paste0("chr",CHR,":",POS,":",OA,":",EA)]

dumIDs = unique(c(myTab_complete$dumID1,myTab_complete$dumID2))
table(is.element(dumIDs,pvar$ID))
dumIDs2 = dumIDs[is.element(dumIDs,pvar$ID)]

myTab_filtered = copy(myTab_complete)
myTab_filtered = myTab_filtered[dumID1 %in% dumIDs2 | dumID2 %in%dumIDs2,]
save(myTab_filtered,file ="../temp/04_CORT_filtered_LIFE.RData")

#' # Best SNP per setting ####
#' ***
myTab_filtered[,absZ := abs(zscore)]

myTab_bestPerSet = copy(myTab_filtered)
myTab_bestPerSet[,dumID3 := paste(setting,cyto,sep = "_")]
setorder(myTab_bestPerSet,-absZ)
myTab_bestPerSet = myTab_bestPerSet[!duplicated(dumID3),]

#' # Save data ####
#' ***
save(myTab_bestPerSet,file ="../results/04_CORT_bestSNPperSetting.RData")

table(is.element(myTab_bestPerSet$dumID2,pvar$ID))

write.table(myTab_bestPerSet[,c(18,6,9)],file = "../temp/04_CORT_ScoreWeights_Combined_best.txt",
            col.names = T,row.names = F,quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

