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

source("../../SourceFile_angmar.R")
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
myTab[SIC %in% D00403[D00403_QUALI_OK==1,SIC], GE := T]
myTab[,table(GE)]
D00403 = D00403[D00403_QUALI_OK==1,]
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
myTab[GENDER==2 & is.na(group) & AGE<50,group := "premenopausal"]

#' ## Sample exclusion criteria
#' 
#' - medication (G03, H02AB, antibaby pill, hormone replacement therapy HRT)
#' - no blood aliquot available
#' - sex in gender questionnaire does not match sex in D00153
#' - women without uterus or ovaries (only relevant in pre-menopausal women)
#' - samples with missing data on time of blood draw, smoking or BMI
#' 
myTab[,goodSample := T]
myTab[D038_G03==T | D038_H02AB==T | D133_med_antibaby==1 | D133_med_HRT==1, goodSample := F]
myTab[is.na(ALIQUOT) & is.na(ALIQUOT_genetics) & is.na(ALIQUOT_GE), goodSample := F]
myTab[GENDER != D133_sex & !is.na(D133_sex), goodSample := F]
myTab[group=="premenopausal" & (D133_removedUterus==1 | D133_removedOvaries==1), goodSample := F]
myTab[is.na(D074_BMI), goodSample := F]
myTab[is.na(D141_smokeStatus), goodSample := F]
myTab[is.na(D126_time), goodSample := F]
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
#' 
myAliquots = unique(c(myTab$ALIQUOT,myTab$ALIQUOT_GE2,myTab$ALIQUOT_genetics))
myAliquots = myAliquots[!is.na(myAliquots)]
myTab2 = data.table(aliquot = myAliquots,
                    GE = F,
                    genetics = F)
myTab2[aliquot %in% myTab$ALIQUOT_GE2,GE := T]
myTab2[aliquot %in% myTab$ALIQUOT_genetics,genetics := T]
matched = match(myTab2$aliquot,myTab$ALIQUOT)
table(is.na(matched))
myTab2[,SIC := myTab[matched,SIC]]
matched = match(myTab2$aliquot,myTab$ALIQUOT_GE2)
table(is.na(matched))
myTab2[,SIC2 := myTab[matched,SIC]]
matched = match(myTab2$aliquot,myTab$ALIQUOT_genetics)
table(is.na(matched))
myTab2[,SIC3 := myTab[matched,SIC]]
myTab2[SIC != SIC2,]
myTab2[SIC != SIC3,]
myTab2[is.na(SIC),SIC:= SIC2]
myTab2[is.na(SIC),SIC:= SIC3]
myTab2[,SIC2 := NULL]
myTab2[,SIC3 := NULL]       
myTab2 = myTab2[,c(4,1:3)]

data2 = data[grepl("T",data)]

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
setnames(myTab2,"ES_V1","E2_S")
setnames(myTab2,"TESTO","TESTO_LCMS")
setnames(myTab2,"PROG","PROG_LCMS")
setnames(myTab2,"ES_MS","E2_LCMS")
setnames(myTab2,"DHEAS","DHEAS_LCMS")
setnames(myTab2,"X17OHP_LCMS","OHP17_LCMS")

save(myTab,myTab2, file = paste0(path_LIFEprepped,"01_LIFEAdult_filtered.RData"))

#' # Check point 2 ####
#' ***
#' Check hormone levels for outlier. I will filter the samples if the value > group mean + 6SD 
#' 
myTab3 = copy(myTab2)
myTab3 = myTab3[!is.na(CORT_S) & !is.na(TESTO_S) & !is.na(E2_S),]
matched = match(myTab3$SIC,myTab$SIC)
table(is.na(matched))
myRows = is.element(names(myTab3),c("aliquot","CORT_S","TESTO_S","E2_S","PROG_LCMS","ALDO_LCMS","DHEAS_S","ANDRO_LCMS","ALDO_LCMS","OHP17_LCMS"))
myRows = seq(1,dim(myTab3)[2],1)[myRows]
myTab3 = cbind(myTab[matched,],myTab3[,myRows,with = F])
names(myTab3)[33] = "ALIQUOT_lab"

myTab3[,QC_ok := T]
myTab3[,reasonToX := ""]

#' ## E2 ####
plot5 = ggplot(myTab3, aes(x=AGE, y=E2_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' I will remove the men with high values (unplausible) - the women with high E2 values might be in the wrong group?
mean = myTab3[,mean(E2_S,na.rm=T),by=group]
sd = myTab3[,sd(E2_S,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "men" & E2_S>400,]
myTab3[group == "men" & E2_S>400,QC_ok := F]
myTab3[group == "men" & E2_S>400,reasonToX := "extreme E2 value"]

myTab3[group == "postmenopausal" & E2_S>556.4,]
myTab3[E2_S>556.4 & group=="postmenopausal",QC_ok := F]
myTab3[E2_S>556.4 & group=="postmenopausal",reasonToX:="high E2 values - consider changing group?"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=E2_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' In the premenopausal women, I want to check how the values match the days since last menstruation
#' 
plot6 = ggplot(myTab3[group=="premenopausal"], aes(x=D133_daysLastMenst, y=E2_S)) +
  #facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("days since last menstruation") + ylab("E2 levels") 
plot6

#' ## P4 ###
plot5 = ggplot(myTab3, aes(x=AGE, y=PROG_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' I will remove the men with high values (unplausible) 
mean = myTab3[,mean(PROG_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(PROG_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "men" & PROG_LCMS>1.8,]
myTab3[group == "men" & PROG_LCMS>1.8,QC_ok := F]
myTab3[group == "men" & PROG_LCMS>1.8,reasonToX := "extreme P4 value"]

myTab3[group == "postmenopausal" & PROG_LCMS>8.7,]
myTab3[group == "postmenopausal" & PROG_LCMS>8.7 & QC_ok ==F ,reasonToX := paste(reasonToX,"extreme P4 value",sep=", ")]
myTab3[group == "postmenopausal" & PROG_LCMS>8.7 & QC_ok==T,reasonToX := "extreme P4 value"]
myTab3[group == "postmenopausal" & PROG_LCMS>8.7,QC_ok := F]

myTab3[group == "premenopausal" & PROG_LCMS>112,]

plot6 = ggplot(myTab3[group=="premenopausal"], aes(x=D133_daysLastMenst, y=PROG_LCMS)) +
  #facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("days since last menstruation") + ylab("P4 levels") 
plot6

#' ## TT ####
plot5 = ggplot(myTab3, aes(x=AGE, y=TESTO_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("TT levels") 
plot5

#' I will remove the two men with high values (unplausible, value > mean + 6SD)
mean = myTab3[,mean(TESTO_S,na.rm=T),by=group]
sd = myTab3[,sd(TESTO_S,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "men" & TESTO_S>54,]
myTab3[group == "men" & TESTO_S>54,QC_ok := F]
myTab3[group == "men" & TESTO_S>54,reasonToX := "extreme TT value"]

myTab3[group == "postmenopausal" & TESTO_S>3.6,]
myTab3[group == "postmenopausal" & TESTO_S>3.6,QC_ok := F]
myTab3[group == "postmenopausal" & TESTO_S>3.6,reasonToX := "extreme TT value"]

myTab3[group == "premenopausal" & TESTO_S>3.7,]

plot5 = ggplot(myTab3[QC_ok==T,], aes(x=AGE, y=TESTO_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("TT levels") 
plot5

#' ## CORT ####
plot5 = ggplot(myTab3, aes(x=AGE, y=CORT_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Cortisol levels") 
plot5

mean = myTab3[,mean(CORT_S,na.rm=T),by=group]
sd = myTab3[,sd(CORT_S,na.rm=T),by=group]
mean$V1 + 6*sd$V1

plot5 = ggplot(myTab3, aes(x=as.factor(D126_time), y=CORT_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_boxplot() +
  theme_bw(base_size = 15) + 
  xlab("Time of blood collection") + ylab("Cortisol levels") 
plot5

#' I want to restrict the analysis to time between 7:00am and 10:59am
myTab3[D126_time==11,]
myTab3[D126_time==11,QC_ok:=F]
myTab3[D126_time==11,reasonToX := "time of blood sampling too late"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=as.factor(D126_time), y=CORT_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_boxplot() +
  theme_bw(base_size = 15) + 
  xlab("Time of blood collection") + ylab("Cortisol levels") 
plot5

#' ## Androstenedione ####
plot5 = ggplot(myTab3, aes(x=AGE, y=ANDRO_LCMS)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Androstenedione levels") 
plot5

#' I will remove the three men & 1 women with high values (unplausible, value > mean + 6SD)
mean = myTab3[,mean(ANDRO_LCMS,na.rm=T),by=GENDER]
sd = myTab3[,sd(ANDRO_LCMS,na.rm=T),by=GENDER]
mean$V1 + 6*sd$V1

myTab3[ANDRO_LCMS>10.7 & GENDER == 1,]
myTab3[ANDRO_LCMS>10.7 & GENDER == 1,QC_ok:=F]
myTab3[ANDRO_LCMS>10.7 & GENDER == 1,reasonToX := "high ANDRO values"]

myTab3[ANDRO_LCMS>8.9 & GENDER == 2,]
myTab3[ANDRO_LCMS>8.9 & GENDER == 2,QC_ok:=F]
myTab3[ANDRO_LCMS>8.9 & GENDER == 2,reasonToX := "high ANDRO values"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=ANDRO_LCMS)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Androstenedione levels") 
plot5

#' ## 17-OHP ####
plot5 = ggplot(myTab3, aes(x=AGE, y=OHP17_LCMS)) +
  facet_wrap(~ GENDER,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("17-OHP levels") 
plot5

mean = myTab3[,mean(OHP17_LCMS,na.rm=T),by=GENDER]
sd = myTab3[,sd(OHP17_LCMS,na.rm=T),by=GENDER]
mean$V1 + 6*sd$V1

myTab3[OHP17_LCMS>17.2 & GENDER == 1,]
myTab3[OHP17_LCMS>17.2 & GENDER == 1 & QC_ok==F,reasonToX := paste(reasonToX,"high 17-OHP values",sep=", ")]
myTab3[OHP17_LCMS>17.2 & GENDER == 1 & QC_ok==T,reasonToX := "high 17-OHP values"]
myTab3[OHP17_LCMS>17.2 & GENDER == 1,QC_ok:=F]

myTab3[OHP17_LCMS>8.1 & GENDER == 2,]
myTab3[OHP17_LCMS>8.1 & GENDER == 2,QC_ok:=F]
myTab3[OHP17_LCMS>8.1 & GENDER == 2,reasonToX := "high 17-OHP values"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=OHP17_LCMS)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("17-OHP levels") 
plot5

#' ## ALDO ###
plot5 = ggplot(myTab3, aes(x=AGE, y=ALDO_LCMS)) +
  facet_wrap(~ GENDER,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Aldosterone levels") 
plot5

mean = myTab3[,mean(ALDO_LCMS,na.rm=T),by=GENDER]
sd = myTab3[,sd(ALDO_LCMS,na.rm=T),by=GENDER]
mean$V1 + 6*sd$V1

myTab3[ALDO_LCMS>715 & GENDER == 1,]
myTab3[ALDO_LCMS>715 & GENDER == 1,QC_ok:=F]
myTab3[ALDO_LCMS>715 & GENDER == 1,reasonToX := "high ALDO values"]

myTab3[ALDO_LCMS>816 & GENDER == 2,]
myTab3[ALDO_LCMS>816 & GENDER == 2,QC_ok:=F]
myTab3[ALDO_LCMS>816 & GENDER == 2,reasonToX := "high ALDO values"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=ALDO_LCMS)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("ALDO levels") 
plot5

#' ## DHEAS ###
plot5 = ggplot(myTab3, aes(x=AGE, y=DHEAS_S)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("DHEA-S levels") 
plot5

mean = myTab3[,mean(DHEAS_S,na.rm=T),by=group]
sd = myTab3[,sd(DHEAS_S,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[DHEAS_S>13.5 & group == "postmenopausal",]
myTab3[DHEAS_S>13.5 & group == "postmenopausal",QC_ok:=F]
myTab3[DHEAS_S>13.5 & group == "postmenopausal",reasonToX := "high DHEAS values"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=DHEAS_S)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("DHEA-S levels") 
plot5


#' # Summary ####
#' ***
myTab3[, table(reasonToX,group)]

#' What is my sample size? 
myTab3[,table(GE, ALIQUOT_GE2==ALIQUOT_lab)]
myTab3[,GE2 := F]
myTab3[ALIQUOT_GE2==ALIQUOT_lab & GE==T,GE2 := T]

myTab3[,table(GE2)]
myTab3[,table(genetics)]
myTab3[,table(genetics,GE2)]

myTab[, table(GE,genetics)]

#' 
#' There are 5168 samples with hormone data (**CORT**, **TT**, and **E2**)
#' 
#' - **TWAS**: There are 2430 samples with hormone AND GE data
#' - **PGS**:  There are 5039 samples with hormone AND genetic data (previous data not stratified for pre- and postmenopausal women)
#' - **TSLS**: There are 2304 samples with hormone AND genetic AND GE data
#' - **eQTL**: There are 2550 samples with genetic AND GE data (previous data not stratified for pre- and postmenopausal women)
#' 
myTab3[, TWAS := GE2]
myTab3[, PGS := genetics]
myTab3[, TSLS := GE2 & genetics]

matched = match(myTab$SIC,myTab3$SIC)
table(is.na(matched))

myTab[,CORT := myTab3[matched,CORT_S]]
myTab[,TESTO := myTab3[matched,TESTO_S]]
myTab[,E2 := myTab3[matched,E2_S]]
myTab[,DHEAS := myTab3[matched,DHEAS_S]]
myTab[,PROG := myTab3[matched,PROG_LCMS]]
myTab[,OHP17 := myTab3[matched,OHP17_LCMS]]
myTab[,ANDRO := myTab3[matched,ANDRO_LCMS]]
myTab[,ALDO := myTab3[matched,ALDO_LCMS]]
myTab[,TWAS := myTab3[matched,TWAS]]
myTab[,PGS := myTab3[matched,PGS]]
myTab[,TSLS := myTab3[matched,TSLS]]
myTab[,ALIQUOT_SH := myTab3[matched,ALIQUOT_lab]]
myTab[,QC_ok := myTab3[matched,QC_ok]]
myTab[,reasonToX := myTab3[matched,reasonToX]]

myTab[,eQTL := GE & genetics]
myTab[,OUremoved := 0]
myTab[,OUremoved1 := myTab3[matched,D133_removedUterus]]
myTab[,OUremoved2 := myTab3[matched,D133_removedOvaries]]
myTab[(OUremoved1==1 & !is.na(OUremoved1)) | (OUremoved2==1 & !is.na(OUremoved2)), OUremoved := 1]
myTab[,OUremoved1 := NULL]
myTab[,OUremoved2 := NULL]
myTab[GENDER==1,OUremoved := NA]

#' Check if there are samples, which are not selected in any setting
#' 
myTab[is.na(eQTL),eQTL := F]
myTab[is.na(TWAS),TWAS := F]
myTab[is.na(PGS),PGS := F]
myTab[is.na(eQTL),eQTL := F]

myTab = myTab[TWAS==T | eQTL==T | TSLS==T | PGS==T,]

#' Filter some columns, which I do not need anymore
#' 
names(myTab)
myNames = c("SIC","EDAT","GRUPPE","group",
            "ALIQUOT","ALIQUOT_genetics","ALIQUOT_GE","ALIQUOT_GE2","ALIQUOT_SH",
            "GENDER","AGE","D126_time","D126_fasting","D074_BMI","D141_smokeStatus",
            "D133_daysLastMenst","OUremoved","CORT","TESTO","E2","DHEAS","PROG","OHP17","ANDRO","ALDO",
            "TWAS","PGS","TSLS","eQTL","QC_ok","reasonToX")
colsOut = setdiff(colnames(myTab),myNames)
myTab[,get("colsOut"):=NULL]
setcolorder(myTab,myNames)
dim(myTab)

#' # Save data ####
#' ***
save(myTab, file = paste0(path_LIFEprepped,"01_LIFEAdult_filtered_final.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

