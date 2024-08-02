#' ---
#' title: "GWAS LIFE"
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
#' Here, I want first to load the pvar files from LIFE and retrict to overlapping SNPs within the predefined candidate gene regions. 
#' 
#' Then I create the necessary phenotype and covariate files to run a GWAS on the selected SNPs and all 8 hormones. 
#'  
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_angmar.R")
.libPaths()

#' # Load genetic data ####
#' ***
load("../results/01_PathwayGenes_byCytoband.RData")
snpdata_Adult = fread(paste0(path_LIFEAdult_Genetics,".snpinfo.txt.gz"))
dim(snpdata_Adult)
snpdata_Heart = fread(paste0(path_LIFEHeart_Genetics,".snpinfo.txt.gz"))
dim(snpdata_Heart)

#' Get overlapping set of SNPs 
snpdata_Adult = snpdata_Adult[id %in% snpdata_Heart$id]
snpdata_Heart = snpdata_Heart[id %in% snpdata_Adult$id]
stopifnot(snpdata_Adult$id == snpdata_Heart$id)
stopifnot(snpdata_Adult$rs == snpdata_Heart$rs)
stopifnot(snpdata_Adult$chr == snpdata_Heart$chr)
stopifnot(snpdata_Adult$position == snpdata_Heart$position)
stopifnot(snpdata_Adult$REF0 == snpdata_Heart$REF0)
stopifnot(snpdata_Adult$ALT1 == snpdata_Heart$ALT1)

#' combine data 
snpdata = copy(snpdata_Adult)
setnames(snpdata,"info","info_A1")
setnames(snpdata,"freq_ALT_imputed","EAF_A1")
setnames(snpdata,"maf_imputed","MAF_A1")
snpdata[,count_ALT_imputed := NULL]
snpdata[,type := NULL]
snpdata[,info_B3 := snpdata_Heart[,info]]
snpdata[,EAF_B3 := snpdata_Heart[,freq_ALT_imputed]]
snpdata[,MAF_B3 := snpdata_Heart[,maf_imputed]]

#' restrict to candidate gene region
dumTab = foreach(i = 1:dim(PathwayGenes_filt)[1])%do%{
  #i=1
  myGene = PathwayGenes_filt[i,]
  message("Working on cytoband ",myGene$cytoband, ", region ",i," of ",dim(PathwayGenes_filt)[1])
  
  data = copy(snpdata)
  data = data[chr == myGene$CHR,]
  data = data[position >= myGene$region_start_collapsed,]
  data = data[position <= myGene$region_end_collapsed,]
  data[,cyto := myGene$cytoband]
  data[,gene := myGene$genes]
  data
}
snpdata_filt = rbindlist(dumTab) 

#' check for duplicates
table(duplicated(snpdata_filt$id), snpdata_filt$cyto)

#' only affects 04q13.2 and 04q13.3
PathwayGenes_filt[grepl("04q13",cytoband)]
snpdata_filt[grepl("04q13",cyto), gene := c("UGT2As | UGT2Bs | SULT1E1")]
snpdata_filt[grepl("04q13",cyto), cyto := "04q13.2-3"]
snpdata_filt = snpdata_filt[!duplicated(id),]
snpdata_filt[,.N,cyto]
save(snpdata_filt,file="../temp/06_SNPList_LIFE.RData")
write.table(snpdata_filt$id,file = "../temp/06_SNPList_LIFE.txt", 
            col.names = F, row.names = F, quote = F)

#' # Create PLINK calls ####
#' ***
#' ## LIFE Adult
#' Create files for PLINK
load(paste0(path_LIFEprepped,"01_LIFEAdult_filtered_final.RData"))
myTab = myTab[QC_ok==T & PGS==T,]

sampleFile = fread(paste0(path_LIFEAdult_Genetics,".psam"))
matched = match(myTab$ALIQUOT_genetics,sampleFile$IID)
table(is.na(matched))

#' These 6 sample were probably filtered for X only, and are now excluded for the autosomes as well
sampleFile = sampleFile[matched,]
sampleFile = sampleFile[!is.na(IID)]
myTab = myTab[ALIQUOT_genetics %in% sampleFile$IID,]
stopifnot(myTab$ALIQUOT_genetics == sampleFile$IID)
stopifnot(myTab$GENDER == sampleFile$SEX)

#' create the phenotype file with 3 * 8 traits
myPhenFile1 = copy(sampleFile)
myPhenFile1[,group := myTab[,group]]

myPhenFile1[group=="men",CORT_M := myTab[group=="men",log(CORT)]]
myPhenFile1[group=="men",TESTO_M := myTab[group=="men",log(TESTO)]]
myPhenFile1[group=="men",E2_M := myTab[group=="men",log(E2)]]
myPhenFile1[group=="men",DHEAS_M := myTab[group=="men",log(DHEAS)]]
myPhenFile1[group=="men",PROG_M := myTab[group=="men",log(PROG)]]
myPhenFile1[group=="men",OHP17_M := myTab[group=="men",log(OHP17)]]
myPhenFile1[group=="men",ANDRO_M := myTab[group=="men",log(ANDRO)]]
myPhenFile1[group=="men",ALDO_M := myTab[group=="men",log(ALDO)]]

myPhenFile1[group=="postmenopausal",CORT_Wpost := myTab[group=="postmenopausal",log(CORT)]]
myPhenFile1[group=="postmenopausal",TESTO_Wpost := myTab[group=="postmenopausal",log(TESTO)]]
myPhenFile1[group=="postmenopausal",E2_Wpost := myTab[group=="postmenopausal",log(E2)]]
myPhenFile1[group=="postmenopausal",DHEAS_Wpost := myTab[group=="postmenopausal",log(DHEAS)]]
myPhenFile1[group=="premenopausal",CORT_Wpre := myTab[group=="premenopausal",log(CORT)]]
myPhenFile1[group=="premenopausal",TESTO_Wpre := myTab[group=="premenopausal",log(TESTO)]]
myPhenFile1[group=="premenopausal",E2_Wpre := myTab[group=="premenopausal",log(E2)]]
myPhenFile1[group=="premenopausal",DHEAS_Wpre := myTab[group=="premenopausal",log(DHEAS)]]

myPhenFile1[group!="men",PROG_W := myTab[group!="men",log(PROG)]]
myPhenFile1[group!="men",OHP17_W := myTab[group!="men",log(OHP17)]]
myPhenFile1[group!="men",ANDRO_W := myTab[group!="men",log(ANDRO)]]
myPhenFile1[group!="men",ALDO_W := myTab[group!="men",log(ALDO)]]

myPhenFile1[,SEX := NULL]
myPhenFile1[,group := NULL]

write.table(myPhenFile1,file=paste0(path_LIFEprepped,"PhenoFile_LIFE_Adult.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

#' create the covariate file with BMI, age and sex
myCovarFile1 = copy(sampleFile)
myCovarFile1[,AGE := myTab$AGE]
myCovarFile1[,BMI := log(myTab$D074_BMI)]
myCovarFile1[,SEX := NULL]

write.table(myCovarFile1,file=paste0(path_LIFEprepped,"CovarFile_LIFE_Adult.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

#' Create PLINK call 
mycall1 = paste0(path_plink2,
                 " --pfile ",path_LIFEAdult_Genetics,
                 " --glm hide-covar firth-fallback",
                 " cols=chrom,pos,ref,alt,firth,test,nobs,machr2,a1freq,a1freqcc,a1countcc,orbeta,se,ci,tz,p",
                 " --pheno ",path_LIFEprepped,"PhenoFile_LIFE_Adult.txt",
                 " --covar ",path_LIFEprepped,"CovarFile_LIFE_Adult.txt",
                 " --extract ../temp/06_SNPList_LIFE.txt",
                 " --out ../results/PLINK_GLM/LIFE_Adult")
mycall1
system(mycall1)

#' ## LIFE Heart
#' Create files for PLINK
load(paste0(path_LIFEprepped,"02_LIFEHeart_filtered_final.RData"))
myTab = myTab[QC_ok==T & PGS==T,]

sampleFile = fread(paste0(path_LIFEHeart_Genetics,".psam"))
sampleFile[,FID:= `#IID`]
sampleFile[,IID := `#IID`]
sampleFile[,`#IID` := NULL]
matched = match(myTab$ALIQUOT_genetics,sampleFile$IID)
table(is.na(matched))

#' These 2 sample were probably filtered for X only, and are now excluded for the autosomes as well
sampleFile = sampleFile[matched,]
sampleFile = sampleFile[!is.na(IID)]
myTab = myTab[ALIQUOT_genetics %in% sampleFile$IID,]
stopifnot(myTab$ALIQUOT_genetics == sampleFile$IID)
stopifnot(myTab$GENDER == sampleFile$SEX)

#' create the phenotype file with 3 * 8 traits
myPhenFile1 = copy(sampleFile)
myPhenFile1[,group := myTab[,group]]

myPhenFile1[group=="men",CORT_M := myTab[group=="men",log(CORT)]]
myPhenFile1[group=="men",TESTO_M := myTab[group=="men",log(TESTO)]]
myPhenFile1[group=="men",E2_M := myTab[group=="men",log(E2)]]
myPhenFile1[group=="men",DHEAS_M := myTab[group=="men",log(DHEAS)]]
myPhenFile1[group=="men",PROG_M := myTab[group=="men",log(PROG)]]
myPhenFile1[group=="men",OHP17_M := myTab[group=="men",log(OHP17)]]
myPhenFile1[group=="men",ANDRO_M := myTab[group=="men",log(ANDRO)]]
myPhenFile1[group=="men",ALDO_M := myTab[group=="men",log(ALDO)]]

myPhenFile1[group=="postmenopausal",CORT_Wpost := myTab[group=="postmenopausal",log(CORT)]]
myPhenFile1[group=="postmenopausal",TESTO_Wpost := myTab[group=="postmenopausal",log(TESTO)]]
myPhenFile1[group=="postmenopausal",E2_Wpost := myTab[group=="postmenopausal",log(E2)]]
myPhenFile1[group=="postmenopausal",DHEAS_Wpost := myTab[group=="postmenopausal",log(DHEAS)]]
myPhenFile1[group=="premenopausal",CORT_Wpre := myTab[group=="premenopausal",log(CORT)]]
myPhenFile1[group=="premenopausal",TESTO_Wpre := myTab[group=="premenopausal",log(TESTO)]]
myPhenFile1[group=="premenopausal",E2_Wpre := myTab[group=="premenopausal",log(E2)]]
myPhenFile1[group=="premenopausal",DHEAS_Wpre := myTab[group=="premenopausal",log(DHEAS)]]

myPhenFile1[group!="men",PROG_W := myTab[group!="men",log(PROG)]]
myPhenFile1[group!="men",OHP17_W := myTab[group!="men",log(OHP17)]]
myPhenFile1[group!="men",ANDRO_W := myTab[group!="men",log(ANDRO)]]
myPhenFile1[group!="men",ALDO_W := myTab[group!="men",log(ALDO)]]

myPhenFile1[,SEX := NULL]
myPhenFile1[,group := NULL]

write.table(myPhenFile1,file=paste0(path_LIFEprepped,"PhenoFile_LIFE_Heart.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

#' create the covariate file with BMI, age and sex
myCovarFile1 = copy(sampleFile)
myCovarFile1[,AGE := myTab$AGE]
myCovarFile1[,BMI := log(myTab$D157_BMI)]
myCovarFile1[,SEX := NULL]

write.table(myCovarFile1,file=paste0(path_LIFEprepped,"CovarFile_LIFE_Heart.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

#' Create PLINK call 
mycall2 = paste0(path_plink2,
                 " --pfile ",path_LIFEHeart_Genetics,
                 " --glm hide-covar firth-fallback",
                 " cols=chrom,pos,ref,alt,firth,test,nobs,machr2,a1freq,a1freqcc,a1countcc,orbeta,se,ci,tz,p",
                 " --pheno ",path_LIFEprepped,"PhenoFile_LIFE_Heart.txt",
                 " --covar ",path_LIFEprepped,"CovarFile_LIFE_Heart.txt",
                 " --extract ../temp/06_SNPList_LIFE.txt",
                 " --out ../results/PLINK_GLM/LIFE_Heart")
mycall2
system(mycall2)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

