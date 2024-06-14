#' ---
#' title: "LIFE-Adult"
#' subtitle: "MR of SH on GE - data preparation"
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
#' In this script, I want to load all the phenotype data I have for LIFE-Adult and bring it to the format I can use throughout the project. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_forostar.R")
.libPaths()

#' # Load derivates data ####
#' ***
data = list.files(path = path_LIFEAdult, pattern = ".xlsx")

#' Okay, there are 25 excel files. I start with D00153, because all samples must have sex and age. 
#' 
D00153 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00153",data)])))
names(D00153)

myTab = copy(D00153)
names(myTab) = gsub("ADULT_PROB_","",names(myTab))
myTab[,dumID := paste(SIC,GRUPPE,sep="_")]
table(duplicated(myTab$SIC))
table(duplicated(myTab$dumID))

#' Now I can add the aliquot ID from the blood sample
D00126 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00126",data)])))
names(D00126)
D00126[,dumID := paste(SIC,GRUPPE,sep="_")]
table(duplicated(D00126$SIC))
table(duplicated(D00126$dumID))
matched = match(myTab$dumID,D00126$dumID)
table(is.na(matched))
D00126 = D00126[matched,]
table(D00126$dumID == myTab$dumID)
table(D00126$ENTNAHME_EDAT == myTab$EDAT)
myTab[,ALIQUOT := D00126$ENTNAHME_SAMPLING_ID]
myTab[,EDAT2 := D00126$ENTNAHME_EDAT]
myTab[,D126_time := D00126[,ENTNAHME_ZP0_START]]
myTab[,D126_time := gsub(":.*","",D126_time)]
myTab[,D126_time := as.numeric(D126_time)]
myTab[,D126_fasting := D00126[,ENTNAHME_FASTED_HOURS]]

#' Okay, now I add a flag for genetic or gene expression data available
#' 
D00364 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00364",data)])))
myTab[,genetics := F]
myTab[SIC %in% D00364$SIC, genetics := T]
myTab[,table(genetics)]
matched = match(myTab$SIC,D00364$SIC)
table(is.na(matched))
myTab[,ALIQUOT_genetics := D00364[matched,ADULT_SNP_SAMPLING_ID]]

D00403 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00403",data)])))
myTab[,GE := F]
myTab[SIC %in% D00403$SIC, GE := T]
myTab[,table(GE)]
matched = match(myTab$SIC,D00403$SIC)
table(is.na(matched))
myTab[,ALIQUOT_GE := D00403[matched,D00403_ALIQUOT]]
myTab[,ALIQUOT_GE2 := gsub("-TR-01","",ALIQUOT_GE)]
myTab[,table(ALIQUOT==ALIQUOT_GE2)]
myTab[,table(ALIQUOT_genetics==ALIQUOT_GE2)]

#' Okay, now I simply match all the other variables (non-blood by SIC, blood parameters by aliquot, must match G)
D00038 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00038",data)])))
D00038[,G03 := grepl("#G03",ADULT_MEDA_H_ATC)]
D00038[,H02AB := grepl("#H02AB",ADULT_MEDA_H_ATC)]
D00038[,table(G03,H02AB)]
matched = match(myTab$SIC,D00038$SIC)
D00038 = D00038[matched,]
myTab[,D038_G03 := D00038[,G03]]
myTab[,D038_H02AB := D00038[,H02AB]]

D00074 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00074",data)])))
matched = match(myTab$SIC,D00074$SIC)
D00074 = D00074[matched,]
myTab[,D074_BMI := D00074[,BMI_BMI]]
myTab[,D074_WHR := D00074[,BMI_WAIST_HIP_RATIO]]

D00141 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00141",data)])))
matched = match(myTab$SIC,D00141$SIC)
D00141 = D00141[matched,]
myTab[,D141_smokeStatus := D00141[,TOB2_SMOKING_STATUS_CURATED]]

D00188 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00188",data)])))
matched = match(myTab$SIC,D00188$SIC)
D00188 = D00188[matched,]
myTab[,D188_diabetesStatus := D00188[,DERIVAT_DIABETES_DS_COMB1]]
myTab[D188_diabetesStatus== -1, D188_diabetesStatus:= NA]
myTab[D188_diabetesStatus== 1, D188_diabetesStatus:= 0]
myTab[D188_diabetesStatus> 1, D188_diabetesStatus:= 1]

D00133 = data.table(read_excel(paste0(path_LIFEAdult,data[grepl("D00133",data)])))
D00133[ADULT_GENDER_GEN01==2,gender := 1]
D00133[ADULT_GENDER_GEN01==1,gender := 2]
matched = match(myTab$SIC,D00133$SIC)
D00133 = D00133[matched,]
myTab[,D133_sex := D00133[,gender]]
myTab[,D133_ageMenarche := D00133[,ADULT_GENDER_FEM01]]
myTab[,D133_daysLastMenst := D00133[,ADULT_GENDER_FEM02]]
myTab[,D133_monthsLastMenst := D00133[,ADULT_GENDER_FEM03]]
myTab[,D133_yearsLastMenst := D00133[,ADULT_GENDER_FEM04]]
myTab[,D133_removedUterus := D00133[,ADULT_GENDER_FEM11]]
myTab[,D133_removedOvaries := D00133[,ADULT_GENDER_FEM13]]
myTab[,D133_med_antibaby := D00133[,ADULT_GENDER_FEM17]]
myTab[,D133_med_HRT := D00133[,ADULT_GENDER_FEM20]]

#' # Check point 1 ####
#' ***
#' ## Menopausal status
#' 
#' I want to get the time since last menstruation: 
#' 
#' - days last menstruation should be between 0-31, -1 is set to NA
#' - month last menstruation should be between 1-12, -1 and 0 is set to NA
#' - years last menstruation should be between 1-x, -1 and 0 is set to NA
#' 
#' Let's define 
#' 
#' - **premenopausal** as all women who are currently menstruating (days not NA) and younger than 60, 
#' - **postmenopausal** as all women which haven't had their menstruation for more than a year (years not NA) and are older than 50 **OR** are older than 60, if years are NA
#' 
myTab[GENDER==1,group := "men"]

myTab[,table(D133_daysLastMenst)]
myTab[D133_daysLastMenst == -1, D133_daysLastMenst := NA]
myTab[,table(D133_yearsLastMenst)]
myTab[D133_yearsLastMenst == -1, D133_yearsLastMenst := NA]

myTab[,table(is.na(D133_yearsLastMenst),AGE>50)]
myTab[GENDER==2 & !is.na(D133_yearsLastMenst) & is.na(D133_daysLastMenst) & D133_yearsLastMenst>=1 & AGE>50,group := "postmenopausal"]
myTab[GENDER==2 & is.na(D133_yearsLastMenst) & AGE>60,group := "postmenopausal"]

myTab[,table(is.na(D133_daysLastMenst),AGE<60)]
myTab[GENDER==2 & !is.na(D133_daysLastMenst) & AGE<60,group := "premenopausal"]

#' ## Sample exclusion criteria
#' 
#' - medication (G03, H02AB, antibaby pill, hormone replacement therapy HRT)
#' - no blood aliquot available
#' - sex in gender questionnaire does not match sex in D00153
#' - women without uterus or ovaries 
#' 
myTab[,goodSample := T]
myTab[D038_G03==T | D038_H02AB==T | D133_med_antibaby==1 | D133_med_HRT==1, goodSample := F]
myTab[is.na(ALIQUOT) & is.na(ALIQUOT_genetics) & is.na(ALIQUOT_GE), goodSample := F]
myTab[GENDER != D133_sex & !is.na(D133_sex), goodSample := F]
myTab[D133_removedUterus==1 | D133_removedOvaries==1, goodSample := F]
myTab[is.na(group), goodSample := F]
myTab[,table(goodSample,group)]

#' ## Save and filter
#' 
save(myTab,file = paste0(path_LIFEprepped,"01_LIFEAdult_unfiltered.RData"))

myTab = myTab[goodSample == T,]

#' ## Exclude remaining duplicates
myTab[,table(duplicated(SIC))]
duplicateSICs = myTab[duplicated(SIC),SIC]
myTab_dups = myTab[SIC %in% duplicateSICs, ]
myTab_dups = myTab_dups[ALIQUOT==ALIQUOT_GE2,]
myTab_nodups = myTab[!is.element(SIC,duplicateSICs), ]
myTab = rbind(myTab_nodups,myTab_dups)

#' # Load blood parameter
#' ***
#' Maybe I can do this per loop? 

data2 = data[grepl("T",data)]

myAliquots = unique(c(myTab$ALIQUOT,myTab$ALIQUOT_GE2,myTab$ALIQUOT_genetics))
myAliquots = myAliquots[!is.na(myAliquots)]
myTab2 = data.table(aliquot = myAliquots,
                    GE = F,
                    genetics = F)
myTab2[aliquot %in% myTab$ALIQUOT_GE2,GE := T]
myTab2[aliquot %in% myTab$ALIQUOT_genetics,genetics := T]

dumTab = foreach(i=1:length(data2))%do%{
  #i=11
  data3 = data.table(read_excel(paste0(path_LIFEAdult,data2[i])))
  if(grepl("DATUM",names(data3)[2])==T){
    parameter1 = gsub("_DATUM","",names(data3)[2])
  }else if(grepl("EDAT",names(data3)[2])==T){
    parameter1 = gsub("_EDAT","",names(data3)[2])
  }
  names(data3) = gsub(paste0(parameter1,"_"),"",names(data3))
  matched = match(myTab2$aliquot,data3$SAMPLING_ID)
  myTab2[,value := data3[matched,NUM_VALUE]]
  setnames(myTab2,"value",parameter1)
}

myTab2[,table(GE,genetics,!is.na(CORT_S))]
names(myTab2)[7] = "E2_S"
names(myTab2)[11] = "TESTO_LCMS"
names(myTab2)[13] = "PROG_LCMS"
names(myTab2)[14] = "E2_LCMS"
names(myTab2)[15] = "DHEAS_LCMS"

save(myTab,myTab2, file = paste0(path_LIFEprepped,"01_LIFEAdult_filtered.RData"))

#' # Check point 2 ####
#' ***
#' Create sample files for each tested scenario
#' 
#' - SH on GE (Aliqout must match!)
#' - genetic on SH (Aliquot must not match - maximize sample size!)
#' - genetic on GE (Aliquot must not match - maximize sample size!)
#' 
dataset1 = copy(myTab2)
dataset1 = dataset1[GE==T,]
dataset1 = dataset1[!is.na(CORT_S) & !is.na(TESTO_S) & !is.na(E2_S),]
matched = match(dataset1$aliquot,myTab$ALIQUOT_GE2)
table(is.na(matched))
dataset1 = cbind(myTab[matched,],dataset1[,c(5,10,7,13,19,6)])

#' Check hormone levels
plot5 = ggplot(dataset1, aes(x=AGE, y=E2_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' I will remove the two men with high values (unplausible) and the four women with values above 250
filt1 = dataset1$group!="premenopausal" & dataset1$E2_S>250
table(filt1)
dataset1 = dataset1[!filt1,]
plot5 = ggplot(dataset1, aes(x=AGE, y=E2_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' In the premenopausal women, I want to check how the values match the days since last menstruation
#' 
plot6 = ggplot(dataset1[group=="premenopausal"], aes(x=D133_daysLastMenst, y=E2_S)) +
  #facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("days since last menstruation") + ylab("E2 levels") 
plot6

plot6 = ggplot(dataset1[group=="premenopausal"], aes(x=D133_daysLastMenst, y=PROG_LCMS)) +
  #facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("days since last menstruation") + ylab("P4 levels") 
plot6

#' How do the other hormones look like
#' 
plot5 = ggplot(dataset1, aes(x=as.factor(D126_time), y=CORT_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_boxplot() +
  theme_bw(base_size = 15) + 
  xlab("Time of blood collection") + ylab("Cortisol levels") 
plot5

plot5 = ggplot(dataset1, aes(x=D074_BMI, y=TESTO_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("BMI") + ylab("Testosterone levels") 
plot5


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

