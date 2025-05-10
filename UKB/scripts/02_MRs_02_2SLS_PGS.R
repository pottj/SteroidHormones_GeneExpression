#' ---
#' title: "Run 2SLS analyses with PGS"
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
library(foreach)
library(doParallel)
library(ivreg)

#' # Parameter settings ####
#' ***
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"
UKB_proteomics = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics/"

#' # Load data ####
#' ***
#' Load phenotype file 
load(paste0(data_QC,"/02_2SLS_input1.RData"))

#' Load SNP info 
load("../temp/proteomics_instruments.RData")
pvar = fread(paste0(data_QC,'01_Prep_02_UKB_Testo_merged.pvar'))
UKB_AF = fread(paste0(data_QC,'01_Prep_02_UKB_Testo_merged_AF.afreq'))

#' Compare allele frequencies 
stopifnot(erg$rsID == UKB_AF$ID)
plot(erg$EAF,UKB_AF$ALT_FREQS)
abline(0,1)
plot(erg$EAF,1-UKB_AF$ALT_FREQS)
abline(0,1)

#' Okay, we always have the different allele. 
#' 
erg[,BETA_harmonized := BETA * (-1)]

#' Create polygenetic score
geno_mat = as.matrix(myTab[,21:28])
betas = as.vector(erg$BETA_harmonized)
PGS = geno_mat %*% betas
myTab[,PGS := PGS]
myTab = myTab[geneticSex == sex,]

# Some tests
model0 = lm(log(TESTO) ~ age +log(BMI), data = myTab)
model2 = lm(log(TESTO) ~ PGS + age + log(BMI) + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10 + Array, data = myTab)
null.r2 = summary(model0)$r.squared
model.r2 = summary(model2)$r.squared
prs.r2 <- model.r2-null.r2
prs.coef <- summary(model2)$coeff["PGS",]
result = data.table(score = "PGS",
                    hormone = "Testosterone", 
                    setting = "male",
                    R2 = prs.r2,
                    sampleSize = length(model2$residuals),
                    beta = prs.coef[1],
                    SE = prs.coef[2],
                    tval = prs.coef[3],
                    pval = prs.coef[4])
result

plotdata = copy(result)

plotdata[,printP := signif(pval, digits = 3)]
plotdata[,printP2 := sub("e", "*x*10^", printP)]
plotdata[,trait := gsub("_", " - ", score)]

myMax = plotdata[,max(R2)]
ggplot(data = plotdata, 
       aes(x = factor(trait), y = R2)) +
  # Specify that we want to print p-value on top of the bars
  geom_text(
    aes(label = paste(printP2)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, myMax * 1.25)) +
  # Specify the axis labels
  xlab("Horomone - sample setting") +
  ylab(expression(paste("PGS model fit:  ", R ^ 2))) +
  # Draw a bar plot
  geom_bar(aes(fill = -log10(pval)), stat = "identity") +
  # Specify the colors
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)
  ) +
  # Some beautification of the plot
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# save the plot

#' Load protein info 
load(paste0(data_QC,"/01_Prep_03_proteinIDs.RData"))

#' Load Olink data
olink = fread(paste0(UKB_proteomics,"olink_data.txt.gz"))

#' # Filter Olink data ####
#' ***
#' Only samples I have in myTab & only probes I have in myGenes
#' 
olink = olink[eid %in% myTab$ID,]
olink = olink[ins_index == 0,]
olink = olink[protein_id %in% myGenes$UKB_coding,]

#' # Loop per protein ####
#' ***
counter = seq(1,dim(myGenes)[1],50)
registerDoParallel(8)

dumTab1 = foreach(i = 1:dim(myGenes)[1])%dopar%{
  #i=1
  myRow = myGenes[i,]
  if(i %in% counter) message("Working on gene ",myRow$GeneSymbol,
                             " (number ",i," of ",dim(myGenes)[1],")")
  
  # filter Olink data
  data1 = copy(olink)
  data1 = data1[protein_id == myRow$UKB_coding,]
  
  # add to main data table
  matched2 = match(myTab$ID,data1$eid)
  myTab[,protein := data1[matched2,result]]
  
  # 2stage least square regression
  ivmodel = ivreg(protein ~ TESTO | PGS, x=TRUE,data = myTab)
  summary(ivmodel)
  ivmodel2 = ivreg(protein ~ TESTO | rs559555 + rs11563251 + rs7694379 + rs4646450 + rs1832007 + rs28929474 + rs2414095 + rs12150660, x=TRUE,data = myTab)
  summary(ivmodel2)
  
  # simple correlation
  cor = cor.test(myTab$TESTO , myTab$protein)
  
  # simple linear regression
  mod1 = lm(protein ~ TESTO + age, data=myTab)
  summary(mod1)
  mod2 = lm(protein ~ TESTO + age + log(BMI), data=myTab)
  summary(mod2)
  
  # summarize
  res = data.table(exposure = "TESTO",
                   outcome = myRow$GeneSymbol,
                   sampleSize = dim(myTab[!is.na(protein)])[1],
                   beta_2SLS1 = c(summary(ivmodel)$coef[2,1]),
                   SE_2SLS1 = c(summary(ivmodel)$coef[2,2]),
                   tval_2SLS1 = c(summary(ivmodel)$coef[2,3]),
                   pval_2SLS1 = c(summary(ivmodel)$coef[2,4]), 
                   beta_2SLS2 = c(summary(ivmodel2)$coef[2,1]),
                   SE_2SLS2 = c(summary(ivmodel2)$coef[2,2]),
                   tval_2SLS2 = c(summary(ivmodel2)$coef[2,3]),
                   pval_2SLS2 = c(summary(ivmodel2)$coef[2,4]), 
                   cor_estimate = cor$estimate,
                   cor_pvalue = cor$p.value,
                   beta_LM = c(summary(mod1)$coef[2,1]),
                   SE_LM = c(summary(mod1)$coef[2,2]),
                   tval_LM = c(summary(mod1)$coef[2,3]),
                   pval_LM = c(summary(mod1)$coef[2,4]),
                   beta_LM2 = c(summary(mod2)$coef[2,1]),
                   SE_LM2 = c(summary(mod2)$coef[2,2]),
                   tval_LM2 = c(summary(mod2)$coef[2,3]),
                   pval_LM2 = c(summary(mod2)$coef[2,4]))
  res
}
TSLS_Tab = rbindlist(dumTab1)
TSLS_Tab[,table(pval_2SLS1<0.05)]
TSLS_Tab[,table(pval_LM<0.05)]
TSLS_Tab[,table(pval_2SLS2<0.05)]
TSLS_Tab[,table(pval_LM2<0.05)]

save(TSLS_Tab,file = "../results/02_2SLS_PGS_results.RData")

myTab[,protein := NULL]
save(myTab,file = paste0(data_QC,"/02_2SLS_input2.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
