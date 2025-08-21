#' ---
#' title: "Get UKB data - sample selection"
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
#' 
#' ## Phenotype data ####
#' 
#' Hormone data: 
#' 
#' - Testosterone: 30850 - 30856
#' - Estradiol: 30800 - 30806
#' 
#' Covariables: 
#' 
#' - sex: 31, 
#' - date of measurement: 53,
#' - UKB center: 54, 
#' - menopause: 2724,
#' - time since last menstrual period: 3700, 
#' - length of menstrual cycle: 3710, 
#' - menstruating today: 3720
#' - medication with steroid hormones: 20003 (48 columns),
#' - smoking status: 20116, 
#' - ethnic background: 21000, 
#' - BMI: 21001, 
#' - age: 21022, 
#' - genetic sex: 22001,
#' - genetic ethnic grouping: 22006,
#' - genetic kinship to other participants: 22021
#' 

myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

#' Get relevant parameters in the annotation file
#'  
hormones = c(paste0("f.3085",c(0:6),".0.0"),paste0("f.3080",c(0:6),".0.0"))
table(is.element(hormones,myAnnot$colNm))

covars = c("f.eid","f.31.0.0","f.53.0.0","f.54.0.0",
           "f.2724.0.0","f.3700.0.0","f.3710.0.0","f.3720.0.0",
           paste0("f.20003.0.",c(0:47)),"f.20116.0.0",
           "f.21000.0.0","f.21001.0.0","f.21022.0.0",
           "f.22001.0.0","f.22006.0.0","f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(hormones,covars),]
myAnnot[,short := c("ID","sex","date","centre","menopause","timeLastPeriod","lengthCycle","menstruatingToday",
                    paste("meds",1:48,sep="_"),"smoking","ethnic","BMI","age",
                    "sex_gen","ethnic_gen","kinship",
                    paste("E2",c("value","date","aliquot","correction","corReason","missReason","reportability"),sep="_"),
                    paste("TESTO",c("value","date","aliquot","correction","corReason","missReason","reportability"),sep="_"))]
myAnnot[,fullDescription := c("Encoded anonymised participant ID",
                              "Sex (data coding 9)",
                              "Date of attending assessment centre",
                              "UK Biobank assessment centre (data coding 10)",
                              "Had menopause (data coding 100579)",
                              "Time since last menstruation (data coding 100291)",
                              "Length of menstrual cycle (data coding 100582)",
                              "Menstruating today (data coding 100349)",
                              rep("Treatment/medication code (data coding 4)",48),
                              "Smoking status *data coding 90)",
                              "Ethnic background (data coding 1001)",
                              "Body mass index",
                              "Age at recruitment", 
                              "Genetic sex (data coding 9)",
                              "Genetic ethnic grouping (data coding 1002)",
                              "Genetic kinship to other participants (data coding 682)",
                              "Oestradiol","Oestradiol assay date",
                              "Oestradiol aliquot (data coding 1529)",
                              "Oestradiol correction level (data coding 165)",
                              "Oestradiol correction reason (data coding 305)",
                              "Oestradiol missing reason (data coding 782)",
                              "Oestradiol reportability (data coding 4917)",
                              "Testosterone","Testosterone assay date",
                              "Testosterone aliquot (data coding 1529)",
                              "Testosterone correction level (data coding 165)",
                              "Testosterone correction reason (data coding 305)",
                              "Testosterone missing reason (data coding 782)",
                              "Testosterone reportability (data coding 4917)")]

x = myAnnot[,colNR]
myTab = fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab)
names(myTab) = myAnnot$short

#' Check for consent (last update: 18/08/2025)
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20250818.csv",UKB_phenotypes))
table(is.element(myTab$ID,ToExclude$V1))
myTab = myTab[!is.element(ID,ToExclude$V1),]

#' Get medication: I want information on steroid medication (ATC starting with G03). I will use the supplemental data 1 file of Wu et al. (Wu, Y., Byrne, E.M., Zheng, Z. et al. Genome-wide association study of medication-use and associated disease in the UK Biobank. Nat Commun 10, 1891 (2019). https://doi.org/10.1038/s41467-019-09572-5)
#' 
MedCoding = data.table(read_xlsx(UKB_MedicationCoding,sheet=1))
hormoneMedication = MedCoding[grepl("G03",ATC_Code),Coding]

myMeds = names(myTab)[grep("meds",names(myTab))]
myTab[,hormoneMeds := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% hormoneMedication,hormoneMeds := 1]
}
myTab[,get("myMeds"):=NULL]
table(myTab$hormoneMeds)

#' Check sex: I only want samples with genetic sex = data base sex. In addition, I want females coded as 2 and males coded as 1
#' 
myTab[,table(sex,sex_gen)]
myTab = myTab[sex == sex_gen,]
myTab[sex==0,sex := 2]

#' ## Proteomics data ####
#' 
#' - PPP sample information: ukb676343.csv (data field 30900 - 30903)
#' - PPP batch: olink_batch_number.dat (see https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=1839, Resources), 
#' - PPP time between blood sampling and measurement: olink_processing_start_date.dat (see https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=1839, Resources)

olink_samples = fread(paste0(UKB_proteomics,"/ukb678244.csv"))
olink_samples = olink_samples[!is.na(`30900-0.0`),]

myTab = myTab[ID %in% olink_samples$eid,]
olink_samples = olink_samples[eid %in% myTab$ID,]
matched = match(myTab$ID,olink_samples$eid)
myTab[,PPP_NRprot := olink_samples[matched,`30900-0.0`],]
myTab[,PPP_plate := olink_samples[matched,`30901-0.0`],]
myTab[,PPP_well := olink_samples[matched,`30902-0.0`],]

olink_batch = fread(paste0(UKB_proteomics,"/olink_batch_number.dat"))
olink_date = fread(paste0(UKB_proteomics,"/olink_processing_start_date.dat"))
matched = match(olink_date$PlateID,olink_batch$PlateID)
olink_date[,batch := olink_batch[matched,Batch]]

matched = match(myTab$PPP_plate,olink_date$PlateID)
table(is.na(matched))
myTab[,PPP_batch := olink_date[matched,batch]]
myTab[,PPP_Startdate := olink_date[matched,Processing_StartDate]]

#' ## Save data 
#' I want to save the phenotype file before any additional filtering. At this moment, it is only filtered for
#' 
#' - consent still given (from 18/08/2025)
#' - sex and genetic sex are the same 
#' - PPP data available
#' 
save(myTab,file = paste0(UKB_dataQC,"/01_UKB_phenotypes_unfiltered.RData"))
# load(paste0(UKB_dataQC,"/01_UKB_phenotypes_unfiltered.RData"))

myAnnot[9,colNm := "f.20003.0.i"]
myAnnot[9,short := "meds_i"]
myAnnot[9,fullDescription := "Treatment/medication code i (i = 0:47, data coding 4)"]
myAnnot = myAnnot[-c(10:56)]

myAnnot2 = data.table(colNm = c(NA,"f.30900.0.0","f.30901.0.0","f.30902.0.0",NA,NA),
                      colNR = c(NA,NA,NA,NA,NA,NA),
                      short = c("hormoneMeds","PPP_NRprot","PPP_plate","PPP_well","PPP_batch","PPP_Startdate"),
                      fullDescription = c("Medication starting with ATC code G03 (summary over all treatment code columns)",
                                          "Pharma Proteomics Project - Number of proteins measured",
                                          "Pharma Proteomics Project - Plate used for sample run",
                                          "Pharma Proteomics Project - Well used for sample run",
                                          "Pharma Proteomics Project - Batch for plate run",
                                          "Pharma Proteomics Project - Start date for plate run (first date per assay)"))
myAnnot = rbind(myAnnot,myAnnot2)
save(myAnnot,file = paste0(UKB_dataQC,"/01_UKB_phenotypes_annotation.RData"))
# load(paste0(UKB_dataQC,"/01_UKB_phenotypes_annotation.RData"))

#' # Check parameters ####
#' ***
#' ## Check covariables ####
#' I want to update some of the variables: 
#' 
#' - menopause: only 0/1 allowed (had menopause? 0 - no, 1 - yes)
#' - time since last menstruation & length of cycle: remove -1, -3. -6 (do not know & prefer not to answer & irregular cycle);
#'    - cycle length between 13-41 (mean + 2SD using UKB showcase) 
#'    - time since last menstruation between 0 and 38 (mean + 1SD using UKB showcase)
#' - smoking: only 0/1 allowed (0 - never & previous, 1 - current)
#' - ethnic: simplify the coding: 
#'    - 1: White, including British (1001), Irish (1002), Any other white background (1003)
#'    - 2: Mixed, including White and Black Caribbean (2001), White and Black African (2002), White and Asian (2003), Any other mixed background (2004)
#'    - 3: Asian or Asian British, including Indian (3001), Pakistani (3002), Bangladeshi (3003), Any other Asian background (3004)
#'    - 4: Black or Black British, including Caribbean (4001), African (4002), Any other Black background (4003)
#'    - 5: Chinese
#'    - 6: Other ethnic groups
#'    
myTab[,table(menopause)]
myTab[menopause %in% c(-3,2,3),menopause := NA]

myTab[,table(timeLastPeriod)]
myTab[timeLastPeriod<0 | timeLastPeriod>38,timeLastPeriod := NA]

myTab[,table(lengthCycle)]
myTab[lengthCycle<13 | lengthCycle>41,lengthCycle := NA]

myTab[,table(smoking)]
myTab[smoking<0,smoking := NA]
myTab[smoking==1,smoking := 0]
myTab[smoking==2,smoking := 1]

myTab[,table(ethnic)]
myTab[,ethnic := substr(ethnic,1,1)]
myTab[ethnic == "-",ethnic := NA]
myTab[,table(ethnic,ethnic_gen)]

hist(myTab$BMI)
hist(myTab$age)

myTab[,sex_gen := NULL]
myTab[,ethnic_gen := NULL]
myTab[,table(kinship)]

myTab[,E2_timeDif := (E2_date - date)/365.25]
myTab[,TESTO_timeDif := (TESTO_date - date)/365.25]
myTab[,PPP_timeDif := (PPP_Startdate - date)/365.25]
hist(myTab$E2_timeDif)
hist(myTab$TESTO_timeDif)
hist(myTab$PPP_timeDif)

myTab[,PPP_plate := gsub("890000000","plate_",PPP_plate)]
myTab[,PPP_batch := paste0("batch_",PPP_batch)]
myTab[,table(PPP_batch)]

save(myTab,file = paste0(UKB_dataQC,"/01_UKB_phenotypes_modified.RData"))

#' ## Check hormones ####
#' I want histograms & boxplots per 5 years for the groups 
#' 
#' - men 
#' - women
#' - post-menopausal women (had menopause == yes) 
#' - pre-menopausal women (had menopause == no & !is.na(timeLastPeriod)
#' 
data1 = copy(myTab)[sex==1,]
data2 = copy(myTab)[sex==2,]
data3 = copy(myTab)[sex==2 & menopause == 1 & age>55,]
data4 = copy(myTab)[sex==2 & menopause == 0 & !is.na(timeLastPeriod) & age<=55,]

data1[,group := "1 - men"]
data2[,group := "2 - women"]
data3[,group := "2b - postmenopausal women"]
data4[,group := "2a - premenopausal women"]

plotData = rbind(data1,data2,data3,data4)
plotData[age<=45 , age_cat := "39-45"]
plotData[age<=50 & age>45, age_cat := "46-50"]
plotData[age<=55 & age>50, age_cat := "51-55"]
plotData[age<=60 & age>55, age_cat := "56-60"]
plotData[age<=65 & age>60, age_cat := "61-65"]
plotData[age<=70 & age>65, age_cat := "66-70"]
plotData[,table(age_cat)]

#' Histograms
plot1.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=log(TESTO_value), col = group,fill=group)) + 
  facet_wrap(~group, scales = "free") +
  geom_histogram(alpha=0.6, binwidth = 0.1)+ 
  xlab("Testosterone levels (log-transformed)") + 
  theme(legend.position = "none")
plot1.1

plot1.2 = ggplot(plotData[!is.na(TESTO_value)], aes(x=TESTO_value, col = group,fill=group)) + 
  facet_wrap(~group, scales = "free") +
  geom_histogram(alpha=0.6, binwidth = 1) + 
  xlab("Testosterone levels") + 
  theme(legend.position = "none")
plot1.2

plot2.1 = ggplot(plotData[!is.na(E2_value)], aes(x=log(E2_value), col = group,fill=group)) + 
  facet_wrap(~group, scales = "free") +
  geom_histogram(alpha=0.6, binwidth = 0.1)+ 
  xlab("Estradiol levels (log-transformed)") + 
  theme(legend.position = "none")
plot2.1

plot2.2 = ggplot(plotData[!is.na(E2_value)], aes(x=E2_value, col = group,fill=group)) + 
  facet_wrap(~group, scales = "free") +
  geom_histogram(alpha=0.6, binwidth = 25) + 
  xlab("Testosterone levels") + 
  theme(legend.position = "none")
plot2.2

#' Boxplot
plot3.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=age_cat, y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("Age (grouped)") + 
  theme(legend.position = "none")
plot3.1

plot3.2 = ggplot(plotData[!is.na(TESTO_value)], aes(x=age_cat, y=TESTO_value, fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels") + 
  xlab("Age (grouped)") + 
  theme(legend.position = "none")
plot3.2

plot4.1 = ggplot(plotData[!is.na(E2_value)], aes(x=age_cat, y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("Age (grouped)") + 
  theme(legend.position = "none")
plot4.1

plot4.2 = ggplot(plotData[!is.na(E2_value)], aes(x=age_cat, y=E2_value, fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels") + 
  xlab("Age (grouped)") + 
  theme(legend.position = "none")
plot4.2

#' save plots 
png(filename = "../results/01_figures/Histogram_TESTO.png",
     width = 2250, height = 1125, res=200)
plot1.1
dev.off()

png(filename = "../results/01_figures/Histogram_logTESTO.png",
    width = 2250, height = 1125, res=200)
plot1.2
dev.off()

png(filename = "../results/01_figures/Histogram_E2.png",
    width = 2250, height = 1125, res=200)
plot2.1
dev.off()

png(filename = "../results/01_figures/Histogram_logE2.png",
    width = 2250, height = 1125, res=200)
plot2.2
dev.off()

png(filename = "../results/01_figures/Boxplot_TESTO.png",
    width = 2250, height = 1125, res=200)
plot3.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logTESTO.png",
    width = 2250, height = 1125, res=200)
plot3.2
dev.off()

png(filename = "../results/01_figures/Boxplot_E2.png",
    width = 2250, height = 1125, res=200)
plot4.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2.png",
    width = 2250, height = 1125, res=200)
plot4.2
dev.off()

#' ## Check other covariables
plot5.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=as.factor(centre), y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("UK Biobank assessment centre") + 
  theme(legend.position = "none")
plot5.1

plot5.2 = ggplot(plotData[!is.na(E2_value)], aes(x=as.factor(centre), y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("UK Biobank assessment centre") + 
  theme(legend.position = "none")
plot5.2

plot6.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=as.factor(smoking), y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("Current smoking (0 - no, 1 - yes)") + 
  theme(legend.position = "none")
plot6.1

plot6.2 = ggplot(plotData[!is.na(E2_value)], aes(x=as.factor(smoking), y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("Current smoking (0 - no, 1 - yes)") + 
  theme(legend.position = "none")
plot6.2

plot7.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=as.factor(ethnic), y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("Ethnic background (1-white, 2-mixed, 3-asian, 4-black, 5-chinese, 6-other)") + 
  theme(legend.position = "none")
plot7.1

plot7.2 = ggplot(plotData[!is.na(E2_value)], aes(x=as.factor(ethnic), y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("Ethnic background (1-white, 2-mixed, 3-asian, 4-black, 5-chinese, 6-other)") + 
  theme(legend.position = "none")
plot7.2

plot8.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=as.factor(kinship), y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("Genetic kinship (-1-excluded, 0-none, 1-1 relative, 10-10 or more relatives)") + 
  theme(legend.position = "none")
plot8.1

plot8.2 = ggplot(plotData[!is.na(E2_value)], aes(x=as.factor(kinship), y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("Genetic kinship (-1-excluded, 0-none, 1-1 relative, 10-10 or more relatives)") + 
  theme(legend.position = "none")
plot8.2

plot9.1 = ggplot(plotData[!is.na(TESTO_value)], aes(x=as.factor(hormoneMeds), y=log(TESTO_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Testosterone levels (log-transformed)") + 
  xlab("Medication with steroid hormones (0 - no, 1 - yes)") + 
  theme(legend.position = "none")
plot9.1

plot9.2 = ggplot(plotData[!is.na(E2_value)], aes(x=as.factor(hormoneMeds), y=log(E2_value), fill=group)) +
  facet_wrap(~group, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  ylab("Estradiol levels (log-transformed)") + 
  xlab("Medication with steroid hormones (0 - no, 1 - yes)") + 
  theme(legend.position = "none")
plot9.2

#' save plots 
png(filename = "../results/01_figures/Boxplot_logTESTO_centre.png",
    width = 2250, height = 1125, res=200)
plot5.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2_centre.png",
    width = 2250, height = 1125, res=200)
plot5.2
dev.off()

png(filename = "../results/01_figures/Boxplot_logTESTO_smoking.png",
    width = 2250, height = 1125, res=200)
plot6.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2_smoking.png",
    width = 2250, height = 1125, res=200)
plot6.2
dev.off()

png(filename = "../results/01_figures/Boxplot_logTESTO_ethnic.png",
    width = 2250, height = 1125, res=200)
plot7.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2_ethnic.png",
    width = 2250, height = 1125, res=200)
plot7.2
dev.off()

png(filename = "../results/01_figures/Boxplot_logTESTO_kinship.png",
    width = 2250, height = 1125, res=200)
plot8.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2_kinship.png",
    width = 2250, height = 1125, res=200)
plot8.2
dev.off()

png(filename = "../results/01_figures/Boxplot_logTESTO_medication.png",
    width = 2250, height = 1125, res=200)
plot9.1
dev.off()

png(filename = "../results/01_figures/Boxplot_logE2_medication.png",
    width = 2250, height = 1125, res=200)
plot9.2
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
