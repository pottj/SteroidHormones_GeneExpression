#' ---
#' title: "Get UKB data - protein selection"
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

#' # Parameter settings ####
#' ***
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"
UKB_proteomics = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics/"
WilsonSupplement = "~/rds/hpc-work/data/downloadedData/2016_Wilson_AndrogenResponseElements_SupTables.xlsx"
eQTLGenData = "~/rds/hpc-work/data/downloadedData/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz"

#' # Load data ####
#' ***
#' ## Load olink data ####
#' ***
olink_coding = fread(paste0(UKB_proteomics,"coding143.tsv"))
olink_coding[,gene:= gsub(";.*","",meaning)]
olink_coding[,description:= gsub(".*;","",meaning)]

#' ## Load Wilson data ####
#' ***
Wilson = data.table(read_excel(WilsonSupplement,sheet=4))
names(Wilson) = c("CHR","START","STOP","Fimo_pval","ARE_Seq","Tier","GenLoc","GeneSymbol","distance_TSS","regulation","up_down")

Wilson = Wilson[regulation == 1,]
setorder(Wilson,"Fimo_pval")
Wilson = Wilson[!duplicated(GeneSymbol),]

#' ## Load eQTLGen data ####
#' ***
eQTLGen = fread(eQTLGenData)
all_genes = eQTLGen[,.N,by=c("Gene","GeneSymbol")]
all_SNPs = eQTLGen[,.N,by = c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele")]
save(all_SNPs,file = "../temp/eQTLGen_allSNPs.RData")
save(all_genes,file = "../temp/eQTLGen_allgenes.RData")
# load("../temp/eQTLGen_allgenes.RData")
# load("../temp/eQTLGen_allSNPs.RData")

load("../temp/proteomics_instruments.RData")
stopifnot(erg$rsID %in% all_SNPs$SNP)

#' ## Load TWAS data ####
#' ***
#' Please check project SteroidHormones_GeneExpression for details! 
#' 
load("../temp/03_ProbeList5_Adult_any.RData")
ProbeList5 = ProbeList5[phenotype == "TESTO",]
ProbeList5 = ProbeList5[setting == "men",]

#' # Merge data ####
#' ***
#' Use eQTLGen as template, and reduce to overlap with olink data
myGenes = copy(all_genes)
myGenes = myGenes[GeneSymbol %in% olink_coding$gene, ]

#' Check for duplicates and remove them (manual look-up in gene cards)
myGenes[duplicated(GeneSymbol)]
myGenes[GeneSymbol %in% c("RNF31","KLRK1","F11R")]
myGenes = myGenes[-grep("ENSG00000259529",Gene)]
myGenes = myGenes[-grep("ENSG00000270149",Gene)]
myGenes = myGenes[-grep("ENSG00000255819",Gene)]

#' Add information of olink and UKB coding 
matched = match(myGenes$GeneSymbol,olink_coding$gene)
myGenes[,UKB_coding := olink_coding[matched,coding]]
myGenes[,description := olink_coding[matched,description]]

#' Add flag for regulation in Wilson data, and add direction of regulation
myGenes[,flag_regulation_Wilson := F]
myGenes[GeneSymbol %in% Wilson$GeneSymbol,flag_regulation_Wilson := T]
matched = match(myGenes$GeneSymbol,Wilson$GeneSymbol)
myGenes[,Wilson_tier := Wilson[matched,Tier]]
myGenes[,Wilson_direction := Wilson[matched,up_down]]
table(myGenes$flag_regulation_Wilson)
myGenes[,Wilson_tier := gsub("two","2",Wilson_tier)]
myGenes[,Wilson_tier := gsub("three","3",Wilson_tier)]
myGenes[,Wilson_tier := gsub("Four","4",Wilson_tier)]
myGenes[,Wilson_tier := gsub("Five","5",Wilson_tier)]
table(myGenes$Wilson_direction,myGenes$Wilson_tier)

#' Add information about significant association in TWAS of LIFE-Adult
myGenes[,flag_association_TWAS := F]
myGenes[GeneSymbol %in% ProbeList5$symbol_INGENUITY,flag_association_TWAS := T]
myGenes[,table(flag_regulation_Wilson,flag_association_TWAS)]

myGenes[flag_regulation_Wilson==T & flag_association_TWAS==T,table(Wilson_tier)]

#' # Save data ####
#' ***
save(myGenes, file = paste0(data_QC,"/01_Prep_03_proteinIDs.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
