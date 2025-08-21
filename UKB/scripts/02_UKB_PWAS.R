#' ---
#' title: "PWAS in UKB"
#' subtitle: "MR of SH on gene expression and protein levels"
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
#' In this script, I want to load the relevant phenotype data for the UK Biobank (UKB). I want to adjust the proteome-wide association studies (PWASs) similar to Sun et al. in their protein GWAS: "For the discovery cohort, association models included the following covariates: age, age2, sex, age × sex, age2 × sex, batch, UKB centre, UKB genetic array, time between blood sampling and measurement and the first 20 genetic principal components." As this is not a genetic study, I will not adjust for UKB genetic array and the first 20 genetic principal components. However, I will check for hormone medication. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")

#' # UKB Data ####
#' ***
load(paste0(UKB_dataQC,"/01_UKB_phenotypes_modified.RData"))
myTab[,logTESTO := log(TESTO_value)]
myTab[,logE2 := log(E2_value)]
myTab[,age2 := age^2]

olink_data = fread(paste0(UKB_proteomics,"/olink_data.txt.gz"))
olink_data = olink_data[ins_index == 0,]
olink_data = olink_data[eid %in% myTab$ID,]

olink_coding = fread(paste0(UKB_proteomics,"/coding143.tsv"))
olink_coding[,gene:= gsub(";.*","",meaning)]
olink_coding[,description:= gsub(".*;","",meaning)]

#' # Loop per protein ####
#' ***
toDoList = data.table(phenotype = rep(c("logTESTO","logE2","TESTO_value","E2_value"),2),
                      transformation = rep(c("log","log","raw","raw"),2),
                      model = c(rep("BMIadj",4),rep("noBMIadj",4)))

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab1 = foreach(i = 1:dim(olink_coding)[1])%dopar%{
#dumTab1 = foreach(i = 1:10)%do%{
  #i=1
  myRow = olink_coding[i,]
  message("Working on ",myRow$gene, " - number ",i," of ",dim(olink_coding)[1])
  myProtein = copy(olink_data)[protein_id == myRow$coding,]
  matched = match(myTab$ID,myProtein$eid)
  myTab[,protein := myProtein[matched,result]]
  
  dumTab2 = foreach(j = 1:dim(toDoList)[1])%do%{
    #j=1
    myToDo = toDoList[j,]
    myTab[,hormone := get(myToDo$phenotype)]
    
    if(myToDo$model == "BMIadj"){
      mod1 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==1)
      mod2 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2)
      mod3 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2 & menopause == 1 & age>55)
      mod4 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + timeLastPeriod + log(BMI) +
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2 & menopause == 0 & !is.na(timeLastPeriod) & age<=55)
    }else{
      mod1 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + #log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==1)
      mod2 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + #log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2)
      mod3 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + #log(BMI) + 
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2 & menopause == 1 & age>55)
      mod4 = lm(protein ~ hormone + age + age2 + smoking + hormoneMeds + timeLastPeriod + #log(BMI) +
                  PPP_timeDif + PPP_batch + as.factor(ethnic) + TESTO_timeDif,
                data=myTab,subset = sex==2 & menopause == 0 & !is.na(timeLastPeriod) & age<=55)
    }
    res2 = data.table(phenotype = rep(myToDo$phenotype,4),
                      setting = c("men","women","postmenopausal","premenopausal"),
                      transformation = rep(myToDo$transformation,4),
                      model = rep(myToDo$model,4),
                      protein_name = rep(myRow$gene,4),
                      protein_id = rep(myRow$coding),
                      N = c(length(mod1$fitted.values),length(mod2$fitted.values),
                            length(mod3$fitted.values),length(mod4$fitted.values)),
                      beta = c(summary(mod1)$coef[2,1],summary(mod2)$coef[2,1],
                               summary(mod3)$coef[2,1],summary(mod4)$coef[2,1]), 
                      SE = c(summary(mod1)$coef[2,2],summary(mod2)$coef[2,2],
                             summary(mod3)$coef[2,2],summary(mod4)$coef[2,2]), 
                      tval = c(summary(mod1)$coef[2,3],summary(mod2)$coef[2,3],
                               summary(mod3)$coef[2,3],summary(mod4)$coef[2,3]), 
                      pval = c(summary(mod1)$coef[2,4],summary(mod2)$coef[2,4],
                               summary(mod3)$coef[2,4],summary(mod4)$coef[2,4])) 
    
    res2
  }
  res1 = rbindlist(dumTab2)
  res1
}

UKB_PWAS = rbindlist(dumTab1) 
UKB_PWAS[,phenotype := gsub("log","",phenotype)]
UKB_PWAS[,phenotype := gsub("_value","",phenotype)]
UKB_PWAS[,dumID := paste(phenotype,setting,model,transformation,sep=":")]
UKB_PWAS[,table(dumID,pval<0.05)]

save(UKB_PWAS,file = "../results/02_UKB_PWAS.RData")
