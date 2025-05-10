#' ---
#' title: "UKB Testosterone in men and women"
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
#' I want to load testosterone levels, restrict to individiuals with Olink-proteom data, and then test for association between testosterone and protein levels per sex. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

library(data.table)
setDTthreads(1)
library(readxl)
library(foreach)
library(doParallel)

#' # Parameter settings ####
#' ***
UKB_SNP_data = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/"
UKB_phenotypes = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/phenotypes/ukb672224.tab"
UKB_proteomics = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics/"
MedicationCoding = "~/rds/hpc-work/data/downloadedData/Wu_2019_SupplementalData1_modified.xlsx"
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"

#' # UKB Data ####
#' ***
myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

#' get relevant parameters: 
#'  
exposure = c("f.30850.0.0")
table(is.element(exposure,myAnnot$colNm))

covars = c("f.eid", "f.31.0.0", paste0("f.20003.0.",c(0:47)),
           "f.21000.0.0","f.21022.0.0","f.22000.0.0",paste0("f.22009.0.",1:10), "f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(exposure,covars),]

x = myAnnot[,colNR]
myTab = fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab)
names(myTab) = c("ID","sex",paste("medication_init",1:48,sep="_"),
                 "ancestry","age","GenotypeBatch",paste("PC",1:10,sep="_"),
                 "kinship","TESTO")

#' ## Proteomics data ####
#' 
olink_samples = fread(paste0(UKB_proteomics,"/ukb676343.csv"))

#' # Filter samples ####
#' ***
#' - testosterone levels available
#' - proteomics available 
#' - White British men without any kinship
#' - No hormone medication (G03)
#' - Consent still true 
#' 
n1 = dim(myTab)[1]

myTab = myTab[!is.na(TESTO),]
n2 = dim(myTab)[1]

myTab = myTab[ID %in% olink_samples[!is.na(`30900-0.0`),eid]]
n3 = dim(myTab)[1]

myTab = myTab[ancestry == 1001,]
myTab = myTab[kinship == 0,]
n4 = dim(myTab)[1]

#' Get medication: I will exclude all men taking testosterone medication (ATC starting with G03)
codingTable = data.table(read_xlsx(MedicationCoding,sheet=1))
hormoneMedication = codingTable[grepl("G03",ATC_Code),Coding]

myMeds = names(myTab)[grep("medication_init",names(myTab))]
myTab[,hormoneMeds := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% hormoneMedication,hormoneMeds := 1]
}
myTab[,get("myMeds"):=NULL]
table(myTab$hormoneMeds)
myTab = myTab[hormoneMeds == 0,]
n5 = dim(myTab)[1]

#' Check for consent
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20241217.csv",UKB_phenotypes))
table(is.element(myTab$ID,ToExclude$V1))
myTab = myTab[!is.element(ID,ToExclude$V1),]
n6 = dim(myTab)[1]

#' Check for outlier
par(mfrow=c(1,2))
hist(myTab[sex==0,TESTO])
hist(myTab[sex==1,TESTO])

threshold_women = myTab[sex==0,mean(TESTO)] + 6*myTab[sex==0,sd(TESTO)] 
threshold_men = myTab[sex==1,mean(TESTO)] + 6*myTab[sex==1,sd(TESTO)] 

myTab = myTab[(sex==0 & TESTO<threshold_women) | (sex==1 & TESTO<threshold_men), ]

boxplot(myTab[sex==0,log(TESTO)],main= "Boxplot of testosterone levels in women",
        sub="after filtering",ylab="Testosterone levels")
hist(myTab[sex==0,log(TESTO)],main= "Histogram of testosterone levels in women",
     sub="after filtering",xlab="Testosterone levels")
boxplot(myTab[sex==1,log(TESTO)],main= "Boxplot of testosterone levels in men",
        sub="after filtering",ylab="Testosterone levels")
hist(myTab[sex==1,log(TESTO)],main= "Histogram of testosterone levels in men",
     sub="after filtering",xlab="Testosterone levels")
myTab[,summary(TESTO),by=sex]

#' # Save phenotype data ####
#' ***
save(myTab, file = paste0(data_QC,"/05_UKB_bothSexes.RData"))

#' # Load Olink data ####
#' ***
olink = fread(paste0(UKB_proteomics,"olink_data.txt.gz"))
olink = olink[eid %in% myTab$ID,]
olink = olink[ins_index == 0,]
olink_coding = fread(paste0(UKB_proteomics,"coding143.tsv"))
olink_coding[,gene:= gsub(";.*","",meaning)]
olink_coding[,description:= gsub(".*;","",meaning)]

#' # Regression Model ####
#' ***
counter = seq(1,dim(olink_coding)[1],100)
#counter = seq(1,100,5)
registerDoParallel(10)
time1<-Sys.time()
dumTab1 = foreach(i = 1:dim(olink_coding)[1])%dopar%{
  #dumTab1 = foreach(i = 1:100)%do%{
  #i=1
  myRow = olink_coding[i,]
  if(i %in% counter) message("Working on gene ",myRow$gene,
                             " (number ",i," of ",dim(olink_coding)[1],")")
  
  # filter Olink data
  data1 = copy(olink)
  data1 = data1[protein_id == myRow$coding,]
  
  # add to main data table
  matched2 = match(myTab$ID,data1$eid)
  myTab[,protein := data1[matched2,result]]
  
  # linear regression
  mod1 = lm(protein ~ log(TESTO) + age, data = myTab, subset = sex==1)
  mod2 = lm(protein ~ log(TESTO) + age, data = myTab, subset = sex==0)
  mod3 = lm(protein ~ (log(TESTO) + age)*sex,data = myTab)
  
  res = data.table(protein = myRow$gene,
                   sex = c("male","female","combined_main","combined_IA"),
                   sampleSize = c(dim(myTab[!is.na(protein) & sex==1])[1],dim(myTab[!is.na(protein) & sex==0])[1],
                                  rep(dim(myTab[!is.na(protein)])[1],2)),
                   beta = c(summary(mod1)$coef[2,1],summary(mod2)$coef[2,1],
                            summary(mod3)$coef[2,1],summary(mod3)$coef[5,1]),
                   SE = c(summary(mod1)$coef[2,2],summary(mod2)$coef[2,2],
                          summary(mod3)$coef[2,2],summary(mod3)$coef[5,2]),
                   tval = c(summary(mod1)$coef[2,3],summary(mod2)$coef[2,3],
                            summary(mod3)$coef[2,3],summary(mod3)$coef[5,3]),
                   pval = c(summary(mod1)$coef[2,4],summary(mod2)$coef[2,4],
                            summary(mod3)$coef[2,4],summary(mod3)$coef[5,4]),
                   adjr2 = c(summary(mod1)$adj.r.squared,summary(mod2)$adj.r.squared,
                             rep(summary(mod3)$adj.r.squared,2)))
  res
}
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")
LR_Tab = rbindlist(dumTab1)
LR_Tab[,table(pval<0.05,sex)]

LR_Tab[,pval_adjusted := p.adjust(pval,method = "fdr"),by=sex]
LR_Tab[,table(pval_adjusted<0.05,sex)]
LR_Tab[,table(sex)]

save(LR_Tab, file = "../results/05_LinReg_MenWomen.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

