#' ---
#' title: "LIFE-Heart"
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
#' In this script, I want to load all the phenotype data I have for LIFE-Heart and bring it to the format I can use throughout the project. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_angmar.R")
.libPaths()

#' # Load derivates data ####
#' ***
data = list.files(path = path_LIFEHeart, pattern = ".xlsx")

#' Okay, there are 9 excel files. I start with D00242, because all samples must have sex and age. 
#' 
D00242 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("D00242",data)])))
names(D00242)

myTab = copy(D00242)
names(myTab) = gsub("HEART_PROB_","",names(myTab))
myTab[,dumID := paste(SIC,GRUPPE,sep="_")]
table(duplicated(myTab$SIC))
table(duplicated(myTab$dumID))

#' Now I can add the aliquot ID from the blood sample
T00991 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("T00991",data)])))
names(T00991)
T00991[,dumID := paste(SIC,GRUPPE,sep="_")]
table(duplicated(T00991$SIC))
table(duplicated(T00991$dumID))
matched = match(myTab$dumID,T00991$dumID)
table(is.na(matched))
T00991 = T00991[matched,]
table(T00991$dumID == myTab$dumID)
table(T00991$HEART_BE_EDAT == myTab$EDAT)
myTab[,ALIQUOT := T00991$HEART_BE_NR]
myTab[,EDAT2 := T00991$HEART_BE_EDAT]
myTab[,T991_time := T00991[,HEART_BE_ZEIT]]
myTab[,T991_time := gsub(":.*","",T991_time)]
myTab[,T991_time := as.numeric(T991_time)]
myTab[,T991_fasting := T00991[,HEART_BE_KARENZZEIT]]

#' Okay, now I add a flag for genetic or gene expression data available
#' 
T00959 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("T00959",data)])))
myTab[,genetics := F]
myTab[SIC %in% T00959$SIC, genetics := T]
myTab[,table(genetics)]
matched = match(myTab$SIC,T00959$SIC)
table(is.na(matched))
myTab[,ALIQUOT_genetics := T00959[matched,B3_SNP_SAMPLINGID]]

T01081 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("T01081",data)])))
myTab[,GE := F]
myTab[SIC %in% T01081[B3_GENEXPR_QUALI_OK==T,SIC], GE := T]
myTab[,table(GE)]
matched = match(myTab$SIC,T01081$SIC)
table(is.na(matched))
myTab[,ALIQUOT_GE := T01081[matched,B3_GENEXPR_SAMPLING_ID]]
myTab[GE==F,ALIQUOT_GE := NA]
myTab[,table(ALIQUOT_genetics==ALIQUOT_GE)]

#' Okay, now I simply match all the other variables (non-blood by SIC, blood parameters by aliquot, must match G)
D00228 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("D00228",data)])))
D00228[,G03 := grepl("#G03",HEART_MEDA_H_ATC)]
D00228[,H02AB := grepl("#H02AB",HEART_MEDA_H_ATC)]
D00228[,table(G03,H02AB)]
matched = match(myTab$SIC,D00228$SIC)
D00228 = D00228[matched,]
myTab[,D228_G03 := D00228[,G03]]
myTab[,D228_H02AB := D00228[,H02AB]]

D00157 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("D00157",data)])))
matched = match(myTab$SIC,D00157$SIC)
D00157 = D00157[matched,]
myTab[,D157_BMI := D00157[,BIOM_ANT_BMI]]
myTab[,D157_WHR := D00157[,BIOM_ANT_WHR]]

D00222 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("D00222",data)])))
matched = match(myTab$SIC,D00222$SIC)
D00222 = D00222[matched,]
myTab[,D222_smokeStatus := D00222[,H_SMO_STATUS]]

D00212 = data.table(read_excel(paste0(path_LIFEHeart,data[grepl("D00212",data)])))
matched = match(myTab$SIC,D00212$pseudonym)
D00212 = D00212[matched,]
myTab[,D121_ageLastMenst := D00212[,nobcad_gender.d00212_fem15]]
myTab[,D121_removedUterus := D00212[,nobcad_gender.d00212_fem05]]
myTab[,D121_removedOvaries := D00212[,nobcad_gender.d00212_fem07]]
myTab[,D121_med_antibaby := D00212[,nobcad_gender.d00212_fem09]]
myTab[,D121_med_HRT := D00212[,nobcad_gender.d00212_fem11]]
myTab[,D121_EDAT_nobCAD := D00212[,edat]]

#' # Check point 1 ####
#' ***
#' ## Menopausal status
#' 
#' In nobcad we asked for age at last menstruation. This should always be smaller or equal to the actual age. But nobcad was done 2014, some years after the baseline. So if age at last menstruation is below the age at baseline, then the proband was post-menopausal in baseline.
#' 
#' Let's define 
#' 
#' - **postmenopausal** as all women which haven't had their menstruation for more than a year (years not NA) and are older than 50 **OR** are older than 60, if years are NA
#' 
myTab[GENDER==1,group := "men"]

#' define age at nobCAD 
myTab[,D121_ageNobCAD := difftime(D121_EDAT_nobCAD,EDAT, units = "weeks")]
myTab[,D121_ageNobCAD := gsub("Time difference of ","",D121_ageNobCAD)]
myTab[,D121_ageNobCAD := gsub(" weeks","",D121_ageNobCAD)]
myTab[,D121_ageNobCAD := as.numeric(D121_ageNobCAD)]
myTab[,D121_ageNobCAD := D121_ageNobCAD/52]
myTab[,D121_ageNobCAD := D121_ageNobCAD + AGE]
plot(myTab$D121_ageNobCAD,myTab$AGE)
abline(0,1)
plot(myTab$D121_ageNobCAD,myTab$D121_ageLastMenst)
abline(0,1)
plot(myTab$AGE,myTab$D121_ageLastMenst)
abline(0,1)
myTab[,D121_yearsLastMenst := AGE - D121_ageLastMenst]
hist(myTab$D121_yearsLastMenst)

#' - years last menstruation should be between 1-x, -1 and 0 is set to NA
#' 

myTab[D121_yearsLastMenst < 1, D121_yearsLastMenst := NA]

myTab[GENDER==2,table(is.na(D121_yearsLastMenst),AGE>50)]
myTab[GENDER==2 & !is.na(D121_yearsLastMenst) & D121_yearsLastMenst>=1 & AGE>50,group := "postmenopausal"]
myTab[GENDER==2 & is.na(D121_yearsLastMenst) & AGE>60,group := "postmenopausal"]

myTab[GENDER==2 & is.na(group) & AGE<60,group := "premenopausal"]

myTab[,table(is.na(group),GE)]
myTab[,table(is.na(group),genetics)]

#' ## Sample exclusion criteria
#' 
#' - medication (G03, H02AB, antibaby pill, hormone replacement therapy HRT)
#' - no blood aliquot available
#' - women without uterus or ovaries (only relevant in pre-menopausal women)
#' - samples with missing data on time of blood draw, smoking or BMI
#' - cohort AMI or special 
#' 
myTab[,goodSample := T]
myTab[D228_G03==T | D228_H02AB==T | D121_med_antibaby==1 | D121_med_HRT==1, goodSample := F]
myTab[is.na(ALIQUOT) & is.na(ALIQUOT_genetics) & is.na(ALIQUOT_GE), goodSample := F]
myTab[group=="premenopausal" & (D121_removedUterus==1 | D121_removedOvaries==1), goodSample := F]
myTab[is.na(D157_BMI), goodSample := F]
myTab[is.na(D222_smokeStatus), goodSample := F]
myTab[is.na(T991_time), goodSample := F]
myTab[is.na(group), goodSample := F]
myTab[COHORT_GROUP %in% c("B3-AMI_01","B3-SPECIAL_01"), goodSample := F]
myTab[,table(goodSample,group)]

#' ## Save and filter
#' 
save(myTab,file = paste0(path_LIFEprepped,"02_LIFEHeart_unfiltered.RData"))

myTab = myTab[goodSample == T,]

#' ## Exclude remaining duplicates
myTab[,table(duplicated(SIC))]
duplicateSICs = myTab[duplicated(SIC),SIC]
myTab_dups = myTab[SIC %in% duplicateSICs, ]
myTab_dups = myTab_dups[ALIQUOT==ALIQUOT_GE,]
myTab_nodups = myTab[!is.element(SIC,duplicateSICs), ]
myTab = rbind(myTab_nodups,myTab_dups)

#' # Load blood parameter
#' ***
#' Maybe I can do this per loop? 
#' 
myAliquots = unique(c(myTab$ALIQUOT,myTab$ALIQUOT_GE,myTab$ALIQUOT_genetics))
myAliquots = myAliquots[!is.na(myAliquots)]
myTab2 = data.table(aliquot = myAliquots,
                    GE = F,
                    genetics = F)
myTab2[aliquot %in% myTab$ALIQUOT_GE,GE := T]
myTab2[aliquot %in% myTab$ALIQUOT_genetics,genetics := T]
matched = match(myTab2$aliquot,myTab$ALIQUOT)
table(is.na(matched))
myTab2[,SIC := myTab[matched,SIC]]
matched = match(myTab2$aliquot,myTab$ALIQUOT_GE)
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

data2 = data[grepl("T0117",data)]

#' **IMPORTANT**: During the LIFE-Heart hormone measurements, there was an update of the LC-MSMS device. Hence, batch 1-4 have different LODs. To make my live simpler, I just kick them all out. 
#' 
dumTab = foreach(i=1:length(data2))%do%{
  #i=1
  data3 = data.table(read_excel(paste0(path_LIFEHeart,data2[i])))
  if(grepl("DATUM",names(data3)[2])==T){
    parameter1 = gsub("_DATUM","",names(data3)[2])
  }else if(grepl("EDAT",names(data3)[2])==T){
    parameter1 = gsub("_EDAT","",names(data3)[2])
  }
  names(data3) = gsub(paste0(parameter1,"_"),"",names(data3))
  data3 = data3[BATCH>=5,]
  matched = match(myTab2$aliquot,data3$SAMPLING_ID)
  myTab2[,value := data3[matched,NUM_VALUE]]
  myTab2[,value := gsub(",",".",value)]
  myTab2[,value := as.numeric(value)]
  setnames(myTab2,"value",parameter1)
  
}
names(myTab2)
myTab2[,table(GE,genetics,!is.na(CORT_LCMS))]
setnames(myTab2,"TESTO_LSMS","TESTO_LCMS")
setnames(myTab2,"ANDRO_LSMS","ANDRO_LCMS")
setnames(myTab2,"OHP_LSMS","OHP17_LCMS")
setnames(myTab2,"ESTR_LCMS","E2_LCMS")

save(myTab,myTab2, file = paste0(path_LIFEprepped,"02_LIFEHeart_filtered.RData"))

#' # Check point 2 ####
#' ***
#' Check hormone levels for outliers. For these plots, I only 
#' 
myTab3 = copy(myTab2)
myTab3 = myTab3[!is.na(CORT_LCMS) & !is.na(TESTO_LCMS) & !is.na(E2_LCMS),]
matched = match(myTab3$SIC,myTab$SIC)
table(is.na(matched))
myRows = is.element(names(myTab3),c("aliquot","CORT_LCMS","TESTO_LCMS","E2_LCMS","PROG_LCMS","ALDO_LCMS","DHEAS_LCMS","ANDRO_LCMS","ALDO_LCMS","OHP17_LCMS"))
myRows = seq(1,dim(myTab3)[2],1)[myRows]
myTab3 = cbind(myTab[matched,],myTab3[,myRows,with = F])
names(myTab3)[31] = "ALIQUOT_lab"

myTab3[,QC_ok := T]
myTab3[,reasonToX := ""]

#' ## E2 ####
plot5 = ggplot(myTab3, aes(x=AGE, y=E2_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

mean = myTab3[,mean(E2_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(E2_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "postmenopausal" & E2_LCMS>580,]
myTab3[E2_LCMS>580 & group=="postmenopausal",QC_ok := F]
myTab3[E2_LCMS>580 & group=="postmenopausal",reasonToX:="high E2 values - consider changing group?"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=E2_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

#' ## TT ####
plot5 = ggplot(myTab3, aes(x=AGE, y=TESTO_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("TT levels") 
plot5

mean = myTab3[,mean(TESTO_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(TESTO_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "postmenopausal" & TESTO_LCMS>4.1,]
myTab3[TESTO_LCMS>4.1 & group=="postmenopausal",QC_ok := F]
myTab3[TESTO_LCMS>4.1 & group=="postmenopausal",reasonToX:="high TT values"]

plot5 = ggplot(myTab3[QC_ok==T,], aes(x=AGE, y=TESTO_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("TT levels") 
plot5

#' ## CORT ####
plot5 = ggplot(myTab3, aes(x=AGE, y=CORT_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Cortisol levels") 
plot5

mean = myTab3[,mean(CORT_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(CORT_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "postmenopausal" & CORT_LCMS>1322,]
myTab3[group == "premenopausal" & CORT_LCMS>1445,]
myTab3[group == "men" & CORT_LCMS>1184,]
myTab3[CORT_LCMS>1184 & group=="men",QC_ok := F]
myTab3[CORT_LCMS>1184 & group=="men",reasonToX:="high CORT values"]

plot5 = ggplot(myTab3, aes(x=as.factor(T991_time), y=CORT_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_boxplot() +
  theme_bw(base_size = 15) + 
  xlab("Time of blood collection") + ylab("Cortisol levels") 
plot5

myTab3[T991_time==13 | T991_time==6,]
myTab3[T991_time==13 | T991_time==6,QC_ok:=F]
myTab3[T991_time==13 | T991_time==6,reasonToX := "time of blood sampling too late/early"]

#' ## A4 ####
plot5 = ggplot(myTab3, aes(x=AGE, y=ANDRO_LCMS)) +
  facet_wrap(~ GENDER,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("Androstenedione levels") 
plot5

mean = myTab3[,mean(ANDRO_LCMS,na.rm=T),by=GENDER]
sd = myTab3[,sd(ANDRO_LCMS,na.rm=T),by=GENDER]
mean$V1 + 6*sd$V1

myTab3[ANDRO_LCMS>9.6 & GENDER == 1,]
myTab3[ANDRO_LCMS>9.6 & GENDER == 1,QC_ok:=F]
myTab3[ANDRO_LCMS>9.6 & GENDER == 1,reasonToX := "high ANDRO values"]

myTab3[ANDRO_LCMS>8.8 & GENDER == 2,]

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

myTab3[OHP17_LCMS>12.3 & GENDER == 1,]
myTab3[OHP17_LCMS>12.3 & GENDER == 1 & QC_ok==F,reasonToX := paste(reasonToX,"high 17-OHP values",sep=", ")]
myTab3[OHP17_LCMS>12.3 & GENDER == 1 & QC_ok==T,reasonToX := "high 17-OHP values"]
myTab3[OHP17_LCMS>12.3 & GENDER == 1,QC_ok:=F]

myTab3[OHP17_LCMS>5.4 & GENDER == 2,]
myTab3[OHP17_LCMS>5.4 & GENDER == 2,QC_ok:=F]
myTab3[OHP17_LCMS>5.4 & GENDER == 2,reasonToX := "high 17-OHP values"]

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

myTab3[ALDO_LCMS>972 & GENDER == 1,]
myTab3[ALDO_LCMS>972 & GENDER == 1,QC_ok:=F]
myTab3[ALDO_LCMS>972 & GENDER == 1,reasonToX := "high ALDO values"]

myTab3[ALDO_LCMS>795 & GENDER == 2,]
myTab3[ALDO_LCMS>795 & GENDER == 2,QC_ok:=F]
myTab3[ALDO_LCMS>795 & GENDER == 2,reasonToX := "high ALDO values"]

plot5 = ggplot(myTab3[QC_ok==T], aes(x=AGE, y=ALDO_LCMS)) +
  facet_wrap(~ as.factor(GENDER),scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("ALDO levels") 
plot5

#' ## DHEAS ###
plot5 = ggplot(myTab3, aes(x=AGE, y=DHEAS_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("DHEA-S levels") 
plot5

mean = myTab3[,mean(DHEAS_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(DHEAS_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[DHEAS_LCMS>9.3 & group == "postmenopausal",]
myTab3[DHEAS_LCMS>17.1 & group == "men",]
myTab3[DHEAS_LCMS>14.7 & group == "premenopausal",]

#' ## P4 ###
plot5 = ggplot(myTab3, aes(x=AGE, y=PROG_LCMS)) +
  facet_wrap(~ group,scales = "free") +
  geom_point() +
  theme_bw(base_size = 15) + 
  xlab("age") + ylab("E2 levels") 
plot5

mean = myTab3[,mean(PROG_LCMS,na.rm=T),by=group]
sd = myTab3[,sd(PROG_LCMS,na.rm=T),by=group]
mean$V1 + 6*sd$V1

myTab3[group == "men" & PROG_LCMS>1.5,]
myTab3[group == "men" & PROG_LCMS>1.5 & QC_ok==F,reasonToX := paste(reasonToX,"extreme P4 value",sep=", ")]
myTab3[group == "men" & PROG_LCMS>1.5 & QC_ok==T,reasonToX := "extreme P4 value"]
myTab3[group == "men" & PROG_LCMS>1.5,QC_ok := F]

myTab3[group == "postmenopausal" & PROG_LCMS>1.6,]
myTab3[group == "postmenopausal" & PROG_LCMS>1.6 & QC_ok ==F ,reasonToX := paste(reasonToX,"extreme P4 value",sep=", ")]
myTab3[group == "postmenopausal" & PROG_LCMS>1.6 & QC_ok==T,reasonToX := "extreme P4 value"]
myTab3[group == "postmenopausal" & PROG_LCMS>1.6,QC_ok := F]

myTab3[group == "premenopausal" & PROG_LCMS>26,]

#' # Summary ####
#' ***
myTab3[, table(reasonToX,group)]

#' What is my sample size? 
myTab3[,table(GE, ALIQUOT_GE==ALIQUOT_lab)]
myTab3[,GE2 := F]
myTab3[ALIQUOT_GE==ALIQUOT_lab,GE2 := T]

myTab3[,table(GE2)]
myTab3[,table(genetics)]
myTab3[,table(genetics,GE2)]

myTab[, table(GE,genetics)]

#' 
#' There are 1821 samples with hormone data (**CORT**, **TT**, and **E2**)
#' - **TWAS**: There are 1676 samples with hormone AND GE data
#' - **PGS**:  There are 1811 samples with hormone AND genetic data 
#' - **TSLS**: There are 1666 samples with hormone AND genetic AND GE data
#' - **eQTL**: There are 2803 samples with genetic AND GE data 
#' 
myTab3[, TWAS := GE2]
myTab3[, PGS := genetics]
myTab3[, TSLS := GE2 & genetics]

matched = match(myTab$SIC,myTab3$SIC)
table(is.na(matched))

myTab[,CORT := myTab3[matched,CORT_LCMS]]
myTab[,TESTO := myTab3[matched,TESTO_LCMS]]
myTab[,E2 := myTab3[matched,E2_LCMS]]
myTab[,DHEAS := myTab3[matched,DHEAS_LCMS]]
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
myTab[(D121_removedUterus==1 & !is.na(D121_removedUterus)) | (D121_removedOvaries==1 & !is.na(D121_removedOvaries)), OUremoved := 1]
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
            "ALIQUOT","ALIQUOT_genetics","ALIQUOT_GE","ALIQUOT_SH",
            "GENDER","AGE","T991_time","T991_fasting","D157_BMI","D222_smokeStatus",
            "OUremoved","CORT","TESTO","E2","DHEAS","PROG","OHP17","ANDRO","ALDO",
            "TWAS","PGS","TSLS","eQTL","QC_ok","reasonToX")
colsOut = setdiff(colnames(myTab),myNames)
myTab[,get("colsOut"):=NULL]
setcolorder(myTab,myNames)
dim(myTab)

#' # Save data ####
#' ***
save(myTab, file = paste0(path_LIFEprepped,"02_LIFEHeart_filtered_final.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

