#' ---
#' title: "PGS - QC of base data"
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
#' Here, I load my five summary statistics file, bring them to the same format and save the filtered data sets. 
#'  
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP150.GRCh38))

#' # Get reference data set ####
#' ***
#' I want to filter for SNPs in both LIFE studies. Filtering for AF and imputation quality will be done later in PLINK. Here, I just filter for overlapping SNP IDs, and no duplicated chrPos
#' 
#' - REF = other allele
#' - ALT = effect allele
#' 
pvar_Adult = fread(paste0(path_LIFEAdult_Genetics,".pvar"))
setnames(pvar_Adult,"#CHROM","CHR")
pvar_Adult[,chrPos := paste0(CHR,":",POS)]
duplicatedSNPs = pvar_Adult[duplicated(chrPos),chrPos]
pvar_Adult = pvar_Adult[!is.element(chrPos,duplicatedSNPs),]

pvar_Heart = fread(paste0(path_LIFEHeart_Genetics,".pvar"))
setnames(pvar_Heart,"#CHROM","CHR")
pvar_Heart[,chrPos := paste0(CHR,":",POS)]
duplicatedSNPs = pvar_Heart[duplicated(chrPos),chrPos]
pvar_Heart = pvar_Heart[!is.element(chrPos,duplicatedSNPs),]

pvar_Adult = pvar_Adult[chrPos %in% pvar_Heart$chrPos,]
pvar_Heart = pvar_Heart[chrPos %in% pvar_Adult$chrPos,]

stopifnot(pvar_Adult$chrPos == pvar_Heart$chrPos)
table(pvar_Adult$ID == pvar_Heart$ID)
pvar_Adult = pvar_Adult[ID %in% pvar_Heart$ID,]
pvar_Heart = pvar_Heart[ID %in% pvar_Adult$ID,]
stopifnot(pvar_Adult$REF == pvar_Heart$REF)
stopifnot(pvar_Adult$ALT == pvar_Heart$ALT)

pvar_Adult[,AF := gsub(";.*","",INFO)]
pvar_Adult[,AF := gsub("AF=","",AF)]
pvar_Adult[,AF := as.numeric(AF)]

pvar_Heart[,AF := gsub(";.*","",INFO)]
pvar_Heart[,AF := gsub("AF=","",AF)]
pvar_Heart[,AF := as.numeric(AF)]

diff = abs(pvar_Adult$AF - pvar_Heart$AF)
filt = diff>0.1
pvar_Adult = pvar_Adult[!filt,]
pvar_Heart = pvar_Heart[!filt,]

pvar = copy(pvar_Adult)
pvar[,FILTER := NULL]
pvar[,INFO := NULL]
setnames(pvar,"AF","AF_Adult")
pvar[,AF_Heart := pvar_Heart$AF]

#' # Create to do list ####
#' ***
#' 
#' ## Testosterone ####
patterns = gsub(".*/","",path_Ruth_TT_2020)
path_Ruth_TT_2020 = gsub(patterns,"",path_Ruth_TT_2020)
statistics = list.files(path = path_Ruth_TT_2020, pattern = patterns)
statistics2 = unlist(strsplit(statistics,"-"))
statistics2 = statistics2[grepl("GCST",statistics2)]

ToDoList_TT = data.table(filenames = statistics,
                      hormone = "TT",
                      setting = c("women","men"),
                      sampleSize = c(230454,194453),
                      path = path_Ruth_TT_2020)

#' ## Estradiol ####
patterns = gsub(".*/","",path_Schmitz_E2_2021)
path_Schmitz_E2_2021 = gsub(patterns,"",path_Schmitz_E2_2021)
statistics = list.files(path = path_Schmitz_E2_2021, pattern = patterns)
statistics2 = unlist(strsplit(statistics,"-"))
statistics2 = statistics2[grepl("GCST",statistics2)]

ToDoList_E2 = data.table(filenames = statistics,
                      hormone = "E2",
                      setting = c("men","women"),
                      sampleSize = c(147690,163985),
                      path = path_Schmitz_E2_2021)

#' ## Cortisol ####
patterns = gsub(".*/","",path_Chan_CORT_2024)
path_Chan_CORT_2024 = gsub(patterns,"",path_Chan_CORT_2024)
statistics = list.files(path = path_Chan_CORT_2024, pattern = patterns)

ToDoList_CORT = data.table(filenames = statistics,
                      hormone = "CORT",
                      setting = c("combined"),
                      sampleSize = c(32981),
                      path = path_Chan_CORT_2024)

#' ## Merge ####
ToDoList = rbind(ToDoList_TT,ToDoList_E2,ToDoList_CORT)

myNames_old_TTE2 = c("hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele",
                     "sampleSize","standard_error","p_value","hm_beta", "hm_effect_allele_frequency",
                     "MAF","phenotype","setting")
myNames_old_CORT = c("rsid","chr","pos","A1","A2",
                     "N","se","P","beta","Freq1",
                     "MAF","phenotype","setting")

myNames_new = c("rsID","CHR","BP","EA","OA","N","SE","P","BETA","EAF","MAF","phenotype","setting")

#' # Load and filter data ####
#' ***
#' At this point, I only want to load the data, extract/rename the relevant columns, and filter for MAF>1% and info>0.8
#' 

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$hormone," in ",myRow$setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(myRow$path,myRow$filenames)
  erg1 = fread(myfn1)
  
  # rename columns
  if(i<5){
    colsOut<-setdiff(colnames(erg1),myNames_old_TTE2)
    erg1[,get("colsOut"):=NULL]
    filt = is.element(myNames_old_TTE2,names(erg1))
    setcolorder(erg1,myNames_old_TTE2[filt])
    names(erg1) = myNames_new[filt]
    
  }else{
    colsOut<-setdiff(colnames(erg1),myNames_old_CORT)
    erg1[,get("colsOut"):=NULL]
    filt = is.element(myNames_old_CORT,names(erg1))
    setcolorder(erg1,myNames_old_CORT[filt])
    names(erg1) = myNames_new[filt]
    
  }
  
  # add missing columns
  myNames_new[!filt]
  if(i<5) erg1[,N := myRow$sampleSize]
  erg1[,MAF := EAF]
  erg1[EAF>0.5,MAF := 1-EAF]
  erg1[,phenotype := myRow$hormone]
  erg1[,setting := myRow$setting]
  stopifnot(sum(names(erg1)%in%myNames_new)==13)
  
  # restrict to data without missing statistics
  erg1 = erg1[!is.na(BETA),]
  erg1 = erg1[!is.na(SE)]
  
  # restrict to data with EAF>1% and EAF<99%
  erg1 = erg1[EAF>=0.01,]
  erg1 = erg1[EAF<=0.99,]
  
  # exclude chromosome X 
  erg1 = erg1[CHR != "X",]
  
  # restrict to rsID available
  erg1 = erg1[grepl("rs",rsID),]
  
  # exclude triallelic SNPs
  erg1[,chrPos := paste(CHR,BP,sep=":")]
  duplicatedSNPs = erg1[duplicated(chrPos),chrPos]
  erg1 = erg1[!is.element(chrPos,duplicatedSNPs),]
  
  # lift over from hg19 to hg38 for cortisol data
  if(i==5){
    mySNPs = erg1$rsID
    snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
    myLift = snpsById(snps, mySNPs, ifnotfound="drop")
    
    myPos = pos(myLift)
    myIDs = myLift$RefSNP_id
    
    matched = match(mySNPs,myIDs)
    erg1[,pos_b38 := myPos[matched]]
    erg1 = erg1[!is.na(pos_b38)]
    erg1[,BP := NULL]
    setnames(erg1,"pos_b38","BP")
    setcolorder(erg1,myNames_new)
    erg1[,chrPos := paste(CHR,BP,sep=":")]
    duplicatedSNPs = erg1[duplicated(chrPos),chrPos]
    erg1 = erg1[!is.element(chrPos,duplicatedSNPs),]
    
  }
  
  # restrict to SNPs in LIFE with matching alleles
  erg1 = erg1[chrPos %in% pvar$chrPos,]
  erg1[,dumID1 := paste0("chr",CHR,":",BP,":",OA,":",EA)]
  erg1[,dumID2 := paste0("chr",CHR,":",BP,":",EA,":",OA)]
  erg1 = erg1[dumID1 %in% pvar$ID | dumID2 %in% pvar$ID,]
  
  matched = match(erg1$chrPos,pvar$chrPos)
  erg1[,REF := pvar[matched, REF]]
  erg1[,ALT := pvar[matched, ALT]]
  erg1[,AF_Adult := pvar[matched, AF_Adult]]
  erg1[,AF_Heart := pvar[matched, AF_Heart]]
  erg1[,flag := F]
  erg1[REF==EA,flag := T]
  
  erg1[flag==T,EA := ALT]
  erg1[flag==T,OA := REF]
  erg1[flag==T,BETA := BETA * (-1)]
  erg1[flag==T,EAF := 1-EAF]
  
  erg1[,diff := abs(EAF - AF_Adult)]
  erg1 = erg1[diff<0.1]
  #erg1[,hist(diff)]
  
  erg1[,ID := paste0("chr",CHR,":",BP,":",OA,":",EA)]
  stopifnot(sum(erg1$ID %in% pvar$ID)==dim(erg1)[1])
  
  # remove some dummy columns
  dumCols = c("chrPos","dumID1","dumID2","REF","ALT","AF_Adult","AF_Heart","flag","diff")
  erg1[,get("dumCols"):=NULL]
  
  # save as txt file and zip it
  outfile = paste0(path_SumStats_QC,myRow$hormone,"_",myRow$setting,".txt")
  fwrite(erg1, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F, na = NA, dec = ".")
  gzip(outfile,destname = paste0(outfile, ".gz"))
  
  # add file name of QC object to todo list
  pvar[,input := F]
  pvar[ID %in% erg1$ID,input := T]
  setnames(pvar,"input",paste0(myRow$hormone,"_",myRow$setting))
  myRow$outfile = paste0(path_SumStats_QC,myRow$hormone,"_",myRow$setting,".txt")
  myRow$SNPs = dim(erg1)[1]
  myRow
}
ToDoList2 = rbindlist(dumTab)

#' # Check ####
#' ***
#' I dont want to filter any SNP at this point, but I want to check that the "shared SNPs" are used with the same alleles. Idealy, they would also have similar effect allele frequency. 
#' 
#' So I will generate a table with chr:pos, and then merge effect allele, other allele, and EAF per study
#' 
flag = rowSums(x = pvar[,c(9:13)])
hist(flag)
table(flag)
pvar[,flag := flag]
pvar = pvar[flag>0,]

write.table(pvar$ID,file = "../results/01_SNPList_raw.txt", 
            col.names = F, row.names = F, quote = F)
save(pvar, file = "../results/01_LIFE_pvar_filtered.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

