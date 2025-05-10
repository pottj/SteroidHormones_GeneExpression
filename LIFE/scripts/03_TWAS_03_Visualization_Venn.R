#' ---
#' title: "TWAS - Visualizing results (Venn diagrams)"
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
#' I want to visualize the overlap between 
#' 
#' - models (BMI adjustment vs. no adjustment): done per hormone and setting and study (10 2-way Venn-Diagrams)
#' - studies (Adult vs Heart): done per hormone and setting (5 4-way Venn-Diagrams) 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/03_TWAS_01_SummaryStatistics_LIFEAdult_hierFDR.RData")
Adult = copy(myTab)

load("../results/03_TWAS_02_SummaryStatistics_LIFEHeart_hierFDR.RData")
Heart = copy(myTab)

stopifnot(names(Heart) == names(Adult))

myTab = rbind(Adult,Heart)

#' # Venn and Euler diagrams in Adult ####
#' ***
#' In Adult, per hormone, setting, and model
#' 
myTab[,uniqueID := paste0(study,"::",uniqueID)]
myTab[,trait := paste0(study,"_",trait)]
myTab[,trait2 := paste0(study,"_",phenotype,"_",setting)]
myTraits = unique(myTab$trait2)
myTraits = myTraits[grepl("Adult",myTraits)]

dumTab = foreach(i = 1:length(myTraits))%do%{
  #i=1
  data1 = copy(myTab)
  data1 = data1[trait2 == myTraits[i],]
  data1 = data1[hierFDR==T,]
  
  IDs_noBMIadj = data1[model=="noBMIadj",PROBE_ID]
  IDs_BMIadj = data1[model=="BMIadj",PROBE_ID]
  
  qlist = venn2(x1 = IDs_noBMIadj,y1 = IDs_BMIadj,
                mytitle = gsub("_"," - ",myTraits[i]),
                mylabels = c("no BMI \nadjustment","BMI \nadjustment"))  
  
  myDat5<-data1[model=="noBMIadj",.(PROBE_ID,beta)]
  myDat6<-data1[model=="BMIadj",.(PROBE_ID,beta)]
  
  makeEuler(mytab1 = myDat5 , mytab2 = myDat6 , 
            object_colname = "PROBE_ID",
            direction_colname = "beta",
            titletext = paste0("Effect direction in ",gsub("_"," - ",myTraits[i])),
            legendtext = c("no BMI \nadjustment","BMI \nadjustment"))
  
  overlap = qlist$q1
  data1 = data1[PROBE_ID %in% overlap,]
  data1
}
myTab2 = rbindlist(dumTab)

#' **Summary**
#'
#' | hormone | sample set | # probes (no BMI) | # probes (BMI adj.) | # probes in overlap | up  | down | 
#' | ------- | ---------- | ----------------: | ------------------: | ------------------: | ---:| ---: |
#' | TT      | men        | 7657              | 281                 | 280                 | 110 | 170  |
#' | E2      | men        | 10                | 10                  | 5                   | 0   | 5    |
#' | TT      | women      | 3                 | 2                   | 0                   | 1   | 1    |
#' | E2      | women      | 1473              | 77                  | 77                  | 18  | 59   |
#' | CORT    | combined   | 973               | 4421                | 930                 | 415 | 515  |
#' 
#' For the time being, I will focus on the genes in the overlap to avoid model problems. 
#' 
#' Now check if these are also associated in LIFE Heart.
#' 
#' # Replication in Heart ####
#' ***
sigProbes = myTab2$uniqueID
sigProbes = gsub("Adult","Heart",sigProbes)
myTab3 = copy(myTab)
myTab3 = myTab3[uniqueID %in% sigProbes,]

dumTab2 = foreach(i = 1:length(myTraits))%do%{
  #i=1
  data1 = copy(myTab2)
  data1 = data1[trait2 == myTraits[i],]
  
  data2 = copy(myTab3)
  data2 = data2[trait2 == gsub("Adult","Heart",myTraits[i]),]
  data2 = data2[P.Value<0.05,]
  
  IDs_Adult = data1[,PROBE_ID]
  IDs_noBMIadj = data2[model=="noBMIadj",PROBE_ID]
  IDs_BMIadj = data2[model=="BMIadj",PROBE_ID]
  
  dummy = gsub("Adult_","",myTraits[i])
  qlist = venn3(x1 = IDs_noBMIadj,y1 = IDs_BMIadj,z1 = IDs_Adult,
                mytitle = gsub("_"," - ",dummy),
                mylabels = c("no BMI \nadjustment","BMI \nadjustment","Adult"))  
  
  myDat5<-data2[model=="noBMIadj" & PROBE_ID %in% qlist$q1,.(PROBE_ID,beta)]
  myDat6<-data2[model=="BMIadj" & PROBE_ID %in% qlist$q1,.(PROBE_ID,beta)]
  myDat7<-data1[model=="noBMIadj" & PROBE_ID %in% qlist$q1,.(PROBE_ID,beta)]
  
  makeEuler(mytab1 = myDat5 , mytab2 = myDat7 , 
            object_colname = "PROBE_ID",
            direction_colname = "beta",
            titletext = paste0("Effect direction in ",gsub("_"," - ",dummy)),
            legendtext = c("Heart","Adult"))
  
  matched = match(myDat5$PROBE_ID,myDat7$PROBE_ID)
  myDat7 = myDat7[matched,]
  stopifnot(myDat5$PROBE_ID == myDat7$PROBE_ID)
  filt = sign(myDat5$beta) == sign(myDat7$beta)
  table(filt)
  myDat5 = myDat5[filt,]
  
  data1 = data1[PROBE_ID %in% myDat5$PROBE_ID,]
  data2 = data2[PROBE_ID %in% myDat5$PROBE_ID,]
  data = rbind(data1,data2)
  data
}
myTab4 = rbindlist(dumTab2)

#' **Summary**
#'
#' | hormone | sample set | # probes in overlap | up  | down | different direction | 
#' | ------- | ---------- | ------------------: | --: | ---: | -------------------:| 
#' | TT      | men        | 53                  | 10  | 40   | 3                   |
#' | E2      | men        | 0                   |     |      |                     |
#' | TT      | women      | 0                   |     |      |                     |
#' | E2      | women      | 2                   | 1   |      | 1                   |
#' | CORT    | combined   | 565                 | 276 | 275  | 14                  |
#' 
#' # Save candidate probes lists ####
#' ***
#' 
#' 1) Probes that can be replicated with the same effect direction in LIFE-Heart (n=602)
#' 2) Probes that can be replicated regardless of BMI adjustment (n=1,292, always same effect direction)
#' 3) Probes that are associated in LIFE-Adult, model without BMI adjustment (n=10,116)
#' 4) Probes that are associated in LIFE-Adult, model with BMI adjustment (n=4,791)
#' 5) Probes that are associated in LIFE-Adult, irrespective of model (n=14,907)
#' 
ProbeList5 = copy(Adult)
ProbeList5 = ProbeList5[hierFDR==T,]

ProbeList4 = copy(Adult)
ProbeList4 = ProbeList4[hierFDR==T & model=="BMIadj",]

ProbeList3 = copy(Adult)
ProbeList3 = ProbeList3[hierFDR==T & model=="noBMIadj",]

ProbeList2 = copy(myTab)
ProbeList2 = ProbeList2[hierFDR==T & uniqueID %in% myTab2$uniqueID,]

ProbeList1 = copy(myTab)
ProbeList1 = ProbeList1[uniqueID %in% myTab4$uniqueID,]

save(ProbeList1, file = "../results/03_TWAS_03_ProbeList1_Replicated_in_Heart.RData")
save(ProbeList2, file = "../results/03_TWAS_03_ProbeList2_Replicated_in_BMImodels.RData")
save(ProbeList3, file = "../results/03_TWAS_03_ProbeList3_Adult_noBMIadj.RData")
save(ProbeList4, file = "../results/03_TWAS_03_ProbeList4_Adult_BMIadj.RData")
save(ProbeList5, file = "../results/03_TWAS_03_ProbeList5_Adult_any.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
