#' ---
#' title: "Get UKB data - sample selection"
#' subtitle: "MR - Testosterone in men on protein levels (Olink)"
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
#' The screening analysis is simply a 2SLS approach: 
#' 
#' protein ~ testosterone | 9 SNPs within/near genes coding for enzyms of the steroid hormone pathway
#' 
#' I want to restrict the analysis to proteins for which I also have gene expression data (eQTLGen phase I). I think it will be interesting to compare these findings later. 
#' 
#' I will add to my gene list the information if they are regulated by an androgen response element (ARE) in prostate cell lines (Wilson et al., 2016), and if they were associated in my TWAS in LIFE-Adult (GE in whole blood)
#' 
#' If there are any significant results, I will attempt a 2-sample MR using Sun et al. for the proteomics data, and Ruth et al. for the testosterone data (higher power because higher sample sizes?)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

library(data.table)
setDTthreads(1)
library(readxl)
library(ggplot2)

#' # Parameter settings ####
#' ***
UKB_SNP_data = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/"
UKB_phenotypes = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/phenotypes/ukb672224.tab"
UKB_proteomics = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics/"
MedicationCoding = "~/rds/hpc-work/data/downloadedData/Wu_2019_SupplementalData1_modified.xlsx"
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"

#' # UKB Data ####
#' ***
#' I want testosterone data and the proteomics data from the BSU UKB application. 
#' 
#' ## Phenotype data ####
#' 
#' - testo            --> 30850
#' - E2               --> 30800
#' 
#' In addition, I want the following columns: 
#' 
#' - ID                 --> eid
#' - Biological sex     --> 31
#' - Medication         --> 20003
#' - Ancestry           --> 21000
#' - BMI                --> 21001
#' - Age                --> 21022
#' - Genotype batch     --> 22000
#' - Genetic PCs        --> 22009
#' - Kinship            --> 22021

myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

#' get relevant parameters: 
#'  
exposure = c("f.30850.0.0","f.30800.0.0")
table(is.element(exposure,myAnnot$colNm))

covars = c("f.eid", "f.31.0.0", paste0("f.20003.0.",c(0:47)),
           "f.21000.0.0","f.21001.0.0","f.21022.0.0","f.22000.0.0",paste0("f.22009.0.",1:10), "f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(exposure,covars),]

x = myAnnot[,colNR]
myTab = fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab)
names(myTab) = c("ID","sex",paste("medication_init",1:48,sep="_"),
                 "ancestry","BMI","age","GenotypeBatch",paste("PC",1:10,sep="_"),
                 "kinship","E2","TESTO")

#' ## Proteomics data ####
#' 
# olink_samples = fread(paste0(UKB_proteomics,"/ukb676343.csv"))

#' # Filter samples ####
#' ***
#' - men
#' - testosterone levels available
#' - proteomics available 
#' - White British men without any kinship, no hormone medication (G03)
#' - Consent still true 
#' 
n0 = dim(myTab)[1]
# myTab = myTab[sex == 1,]
# n1 = dim(myTab)[1]

# myTab = myTab[!is.na(TESTO) | !is.na(E2),]
# n2 = dim(myTab)[1]
# myTab[!is.na(E2)]
# myTab = myTab[ID %in% olink_samples[!is.na(`30900-0.0`),eid]]
# n3 = dim(myTab)[1]

myTab = myTab[ancestry == 1001,]
myTab = myTab[kinship == 0,]

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
n4 = dim(myTab)[1]

#' Check for consent
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20241217.csv",UKB_phenotypes))
table(is.element(myTab$ID,ToExclude$V1))
myTab = myTab[!is.element(ID,ToExclude$V1),]
n5 = dim(myTab)[1]

#' # Plotting ####
#' ***
myTab6 = copy(myTab)

myTab6[age<=45 , age3 := "38-45"]
myTab6[age<=50 & age>45, age3 := "46-50"]
myTab6[age<=55 & age>50, age3 := "51-55"]
myTab6[age<=60 & age>55, age3 := "56-60"]
myTab6[age<=65 & age>60, age3 := "61-65"]
myTab6[age>65, age3 := "66-72"]
myTab6[,table(age3)]

myTab6[sex==0, sex2 := "female"]
myTab6[sex==1, sex2 := "male"]
myTab6[,table(age3,sex2)]

myTab6 = myTab6[QC_ok==T,]

#' ## Estradiol
ggp2 = ggplot(myTab6, aes(x=age3, y=E2, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Estradiol levels (in pmol/l)", 
       fill="Sex",
       title = paste0("UK Biobank")) +
  theme_classic() +
  theme(legend.position = "none")
ggp2

ggp3 = ggplot(myTab6[age>50], aes(x=age3, y=E2, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Estradiol levels (in pmol/l)", 
       fill="Sex",
       title = paste0("UK Biobank (age>50y)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp3

#' ## Testosterone
ggp4 = ggplot(myTab6, aes(x=age3, y=TESTO, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Testosterone levels (in nmol/l)", 
       fill="Sex",
       title = paste0("UK Biobank")) +
  theme_classic() +
  theme(legend.position = "none")
ggp4

ggp5 = ggplot(myTab6[age>50,], aes(x=age3, y=TESTO, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Testosterone levels (in nmol/l)", 
       fill="Sex",
       title = paste0("UK Biobank (age>50y)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp5

#' ## Testosterone / Estradiol
myTab6[,ratio := TESTO*1000/E2]

ggp6 = ggplot(myTab6, aes(x=age3, y=ratio, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Testosterone to Estradiol ratio", 
       fill="Sex",
       title = paste0("UK Biobank")) +
  theme_classic() +
  theme(legend.position = "none")
ggp6

ggp7 = ggplot(myTab6[age>50,], aes(x=age3, y=ratio, fill=sex2)) +
  facet_wrap(~sex2, scales = "free_y") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Testosterone to Estradiol ratio", 
       fill="Sex",
       title = paste0("UK Biobank (age>50y)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp7

#' # Save plots ####
#' ***
ggsave(plot = ggp2, filename = paste0("../results/_figures/Boxplots_UKB_Estradiol.png"), 
       height = 7, width = 14)

ggsave(plot = ggp3, filename = paste0("../results/_figures/Boxplots_UKB_Estradiol_50y.png"), 
       height = 7, width = 14)

ggsave(plot = ggp4, filename = paste0("../results/_figures/Boxplots_UKB_Testosterone.png"), 
       height = 7, width = 14)

ggsave(plot = ggp5, filename = paste0("../results/_figures/Boxplots_UKB_Testosterone_50y.png"), 
       height = 7, width = 14)

ggsave(plot = ggp6, filename = paste0("../results/_figures/Boxplots_UKB_Testosterone_Estradiol.png"), 
       height = 7, width = 14)

ggsave(plot = ggp7, filename = paste0("../results/_figures/Boxplots_UKB_Testosterone_Estradiol_50y.png"), 
       height = 7, width = 14)
