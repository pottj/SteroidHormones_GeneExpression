#' ---
#' title: "2SLS analyses evaluation"
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
library(biomaRt)

#' # Parameter settings ####
#' ***
data_QC = "~/rds/hpc-work/MR_Testo_Proteomics/"

#' # Load data ####
#' ***
load("../results/02_2SLS_PGS_results.RData")
TSLS_Tab[,pval_2SLS_adj := p.adjust(pval_2SLS,method = "fdr")]
candidateGenes = TSLS_Tab[pval_2SLS_adj <0.05,outcome]

load("../temp/proteomics_instruments.RData")
load(paste0(data_QC,"/01_Prep_03_proteinIDs.RData"))

#' # Add gene information ####
#' ***
ensembl = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
myWishList = c('ensembl_gene_id', "description", 
               'hgnc_symbol','uniprotswissprot',"entrezgene_id",
               "band", 'chromosome_name','start_position','end_position',
               "transcript_start","transcript_end","transcription_start_site",
               "gene_biotype")
all_genes = getBM(attributes=myWishList, mart = ensembl)
head(all_genes)
setDT(all_genes)

candidate_genes = all_genes[ensembl_gene_id %in% myGenes[GeneSymbol %in% candidateGenes,Gene],]
candidate_genes = candidate_genes[uniprotswissprot != ""]
candidate_genes = candidate_genes[!duplicated(ensembl_gene_id),]
candidate_genes[,cytoband := paste0(chromosome_name,band)]
candidate_genes[chromosome_name %in% c(1:9),cytoband := paste0("0",chromosome_name,band)]
candidate_genes[,flag := T]
candidate_genes[cytoband %in% erg$cyto,flag := F]

#' # Combine and save ####
#' ***
topList = TSLS_Tab[pval_2SLS_adj <0.05,]
matched = match(topList$outcome,candidate_genes$hgnc_symbol)
topList[, cytoband := candidate_genes[matched,cytoband]]
topList[, description := candidate_genes[matched,description]]
topList[, uniprot := candidate_genes[matched,uniprotswissprot]]
topList[, start_hg38 := candidate_genes[matched,start_position]]
topList[, end_hg38 := candidate_genes[matched,end_position]]
topList[, pleiotropy := F]
topList[cytoband %in% erg$cyto, pleiotropy := T]
topList[pleiotropy == F, ]

save(topList,file = "../results/03_TopList_2SLS_PGS.RData")

erg[,BP_hg19 := c(31810974, 234679384 , 88186509, 99266318,
                 5254847 , 94844947 , 51524292 , 7521915 )]
save(erg,file = "../temp/proteomics_instruments_hg19.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
