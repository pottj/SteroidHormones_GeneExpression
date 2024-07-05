#' ---
#' title: "Get E2 SNPs"
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
#' In this script, I will load the estradiol data from Schmitz et al.(2021), and filter for 
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

source("../../SourceFile_forostar.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/01_PathwayGenes_byCytoband.RData")

patterns = gsub(".*/","",path_Schmitz_E2_2021)
path_Schmitz_E2_2021 = gsub(patterns,"",path_Schmitz_E2_2021)
statistics = list.files(path = path_Schmitz_E2_2021, pattern = patterns)

statistics2 = unlist(strsplit(statistics,"-"))
statistics2 = statistics2[grepl("GCST",statistics2)]

ToDoList = data.table(filenames = statistics,
                      hormone = "E2",
                      setting = c("men","women"),
                      sampleSize = c(147690,163985))

myNames_old = c("hormone","setting",
                "hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_effect_allele_frequency",
                "hm_beta","standard_error","p_value","zscore","sampleSize")
myNames_new = c("hormone","setting",
                "rsID","CHR","POS","EA","OA","EAF",
                "beta","SE","pvalue","zscore","sampleSize")

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myTime = Sys.time()
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$hormone," in ",myRow$setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(path_Schmitz_E2_2021,myRow$filenames)
  erg1 = fread(myfn1)
  erg1 = erg1[hm_pos == base_pair_location,]
  erg1[is.na(hm_beta),hm_rsid := variant_id]
  erg1[is.na(hm_beta),hm_chrom := chromosome]
  erg1[is.na(hm_beta),hm_pos := base_pair_location]
  erg1[is.na(hm_beta),hm_effect_allele := effect_allele]
  erg1[is.na(hm_beta),hm_other_allele := other_allele]
  erg1[is.na(hm_beta),hm_effect_allele_frequency := effect_allele_frequency]
  erg1[is.na(hm_beta),hm_beta := beta]
  
  erg1[,hormone := myRow$hormone]
  erg1[,setting := myRow$setting]
  erg1[,sampleSize := myRow$sampleSize]
  erg1[,zscore := hm_beta/standard_error]
  
  stopifnot(myNames_old %in% names(erg1))
  colsOut<-setdiff(colnames(erg1),myNames_old)
  erg1[,get("colsOut"):=NULL]
  setcolorder(erg1,myNames_old)
  names(erg1) = myNames_new
  
  # get list 3: all SNPs of the pathway gene regions
  dumTab2 = foreach(j = 1:dim(PathwayGenes_filt)[1])%do%{
    #j=1
    message("Working on locus number ",j)
    erg4 = copy(erg1)
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

save(myTab_complete,file ="../temp/03_E2_unfiltered.RData")

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
save(myTab_filtered,file ="../temp/03_E2_filtered_LIFE.RData")

#' # Best SNP per setting ####
#' ***
myTab_filtered[,absZ := abs(zscore)]

myTab_bestPerSet = copy(myTab_filtered)
myTab_bestPerSet[,dumID3 := paste(setting,cyto,sep = "_")]
setorder(myTab_bestPerSet,-absZ)
myTab_bestPerSet = myTab_bestPerSet[!duplicated(dumID3),]

#' # Same SNP set ####
#' ***
myTab_sameSNPSet = copy(myTab_filtered)
cytoList2 = myTab_sameSNPSet[absZ>5,.N,by = c("cyto","setting")]
cytoList3 = dcast(cytoList2,cyto ~ setting,value.var = "N")
setDT(cytoList3)

dumTab3 = foreach(j=1:dim(cytoList3)[1])%do%{
  #j=1
  myTab_M = copy(myTab_sameSNPSet)[setting == "men" & cyto == cytoList3[j,cyto],]
  myTab_W = copy(myTab_sameSNPSet)[setting == "women" & cyto == cytoList3[j,cyto],]
  
  dups = myTab_W[duplicated(rsID),rsID] 
  myTab_W = myTab_W[!is.element(rsID,dups)]
  dups = myTab_M[duplicated(rsID),rsID] 
  myTab_M = myTab_M[!is.element(rsID,dups)]
  
  myTab_M = myTab_M[rsID %in% myTab_W$rsID]
  myTab_W = myTab_W[rsID %in% myTab_M$rsID]
  
  dataset_women = list(beta=myTab_W[,beta], 
                       varbeta=(myTab_W[,SE])^2,
                       N=myTab_W[,sampleSize],
                       snp=myTab_W[,rsID],
                       type="quant",
                       MAF = myTab_W[,MAF],
                       position = myTab_W[,POS])
  dataset_men = list(beta=myTab_M[,beta], 
                     varbeta=(myTab_M[,SE])^2,
                     N=myTab_M[,sampleSize],
                     snp=myTab_M[,rsID],
                     type="quant",
                     MAF = myTab_M[,MAF],
                     position = myTab_M[,POS])
  
  my_res1 = coloc::coloc.abf(dataset1=dataset_women,dataset2=dataset_men)
  my_res2 = my_res1$summary
  x2<-as.data.table(my_res2)
  x3<-t(x2)
  x4<-as.data.table(x3)
  names(x4)<-names(my_res2)
  x4[,cytoband:=cytoList3[j,cyto]]
  x4[,trait1:="women"]
  x4[,trait2:="men"]
  snp1 = myTab_W[abs(zscore) == max(abs(zscore)),rsID]
  x4[,bestSNP_trait1 := snp1[1]]
  x4[,bestSNP_trait1_zscore := myTab_W[rsID == snp1[1],abs(zscore)]]
  snp2 = myTab_M[abs(zscore) == max(abs(zscore)),rsID]
  x4[,bestSNP_trait2 := snp2[1]]
  x4[,bestSNP_trait2_zscore := myTab_M[rsID == snp2[1],abs(zscore)]]
  x4
  
}
myColoc = rbindlist(dumTab3)
stopifnot(myColoc$cytoband == cytoList3$cyto)
myColoc[PP.H4.abf>0.5,comment := "shared causal SNP"]
myColoc[PP.H3.abf>0.5,comment := "two independent SNPs"]
myColoc[PP.H2.abf>0.5,comment := "only in men"]
myColoc[PP.H1.abf>0.5,comment := "only in women"]
table(myColoc$comment)
myColoc = cbind(myColoc,cytoList3[,c(2:3),with=F])

dumID1 = myColoc[PP.H1.abf>0.5,bestSNP_trait1]
dumID2 = myColoc[PP.H2.abf>0.5,bestSNP_trait2]
dumID31 = myColoc[PP.H3.abf>0.5,bestSNP_trait1]
dumID32 = myColoc[PP.H3.abf>0.5,bestSNP_trait2]
dumID41 = myColoc[PP.H4.abf>0.5 & bestSNP_trait2_zscore<bestSNP_trait1_zscore,bestSNP_trait1]
dumID42 = myColoc[PP.H4.abf>0.5 & bestSNP_trait2_zscore>bestSNP_trait1_zscore,bestSNP_trait2]

dumIDs = c(dumID1,dumID2,dumID31,dumID32,dumID41,dumID42)

#' Now I can filter for those selected instruments
myTab_sameSNPSet = myTab_sameSNPSet[rsID %in% dumIDs,]
myTab_sameSNPSet[rsID %in% dumID1, flag := "only in women"]
myTab_sameSNPSet[rsID %in% dumID2, flag := "only in men"]
myTab_sameSNPSet[rsID %in% dumID31 | rsID %in% dumID32, flag := "two independent SNPs"]
myTab_sameSNPSet[rsID %in% dumID41 | rsID %in% dumID42, flag := "shared causal SNP"]
filt = myTab_sameSNPSet$cyto == "04q13.2" & myTab_sameSNPSet$rsID == "rs4464622"
table(filt)
myTab_sameSNPSet = myTab_sameSNPSet[!filt, ]

#' # Save data ####
#' ***
save(myTab_bestPerSet,file ="../results/03_E2_bestSNPperSetting.RData")
save(myTab_sameSNPSet,file = "../results/03_E2_sharedOrIndepSNPs.RData")
save(myColoc,file = "../results/03_E2_Coloc.RData")

table(is.element(myTab_bestPerSet$dumID2,pvar$ID))
table(is.element(myTab_sameSNPSet$dumID2,pvar$ID))

write.table(myTab_bestPerSet[setting=="women",c(18,6,9)],file = "../temp/03_E2_ScoreWeights_Women_best.txt",
            col.names = T,row.names = F,quote = F)
write.table(myTab_bestPerSet[setting=="men",c(18,6,9)],file = "../temp/03_E2_ScoreWeights_Men_best.txt",
            col.names = T,row.names = F,quote = F)

write.table(myTab_sameSNPSet[setting=="women",c(18,6,9)],file = "../temp/03_E2_ScoreWeights_Women_same.txt",
            col.names = T,row.names = F,quote = F)
write.table(myTab_sameSNPSet[setting=="men",c(18,6,9)],file = "../temp/03_E2_ScoreWeights_Men_same.txt",
            col.names = T,row.names = F,quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

