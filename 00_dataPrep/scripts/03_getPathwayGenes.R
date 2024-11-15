#' ---
#' title: "Steroid hormone pathway genes"
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
#' I downloaded from [KEGG](https://www.genome.jp/dbget-bin/www_bget?path:hsa00140) the HGNC symbol for 62 genes involved in the steroid hormone pathway. This includes the KEGG entry ID, KEGG orthology ID, the gene symbol, and the full gene name. 
#' 
#' I add *SHBG* and *CBG* (also known as *SERPINA6*) manually. They are not part of the pathway but transport the hormones from production site to effect site. 
#' 
#' Here, I use bioMart to add the chromosome, start and end position (GRCh38), and cytoband of the gene. 
#' 
#' Finally, I collapse neighboring genes (most likely when positioned at the same cytoband).
#'  
#' **Note**: depending on admin settings, this script might not run in RStudio as it needs the server to connect to Ensembl. 
#'  
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "angmar"

source("../../SourceFile.R")
.libPaths()

#' # Load gene list ####
#' ***
PathwayGenes = fread(path_PathwayGenesKEGG)
PathwayGenes[1:5,]

SHBG = data.table(KEGG_entryID = 6462,
                  KEGG_orthologyID = "K25754",
                  gene = "SHBG",
                  gene_fullName = "sex hormone-binding globulin")

CBG = data.table(KEGG_entryID = 866,
                 KEGG_orthologyID = "K04525",
                 gene = "SERPINA6",
                 gene_fullName = "serpin family A member 6, CBG")

PathwayGenes = rbind(PathwayGenes,SHBG,CBG,fill=T)

#' # Get gene information ####
#' ***
biomartCacheClear()
ensembl_38 = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")

myWishList = c('ensembl_gene_id', "description", 
               'hgnc_symbol','uniprotswissprot',"entrezgene_id",
               "band", 'chromosome_name','start_position','end_position',
               "transcript_start","transcript_end","transcription_start_site",
               "gene_biotype")

all_genes_b38 = getBM(attributes=myWishList, mart = ensembl_38)
setDT(all_genes_b38)
all_genes_b38 = all_genes_b38[hgnc_symbol %in% PathwayGenes$gene,]
write.table(all_genes_b38, "../results/03_biomart_genes_GRCh38.txt",row.names = F)

all_genes_b38[ ,table(chromosome_name)]
all_genes_b38 = all_genes_b38[chromosome_name %in% c(1:22,"X"),]

test1 = all_genes_b38[,unique(chromosome_name),by=hgnc_symbol]
matched = match(PathwayGenes$gene, test1$hgnc_symbol)
test1 = test1[matched,]

test2 = all_genes_b38[,unique(start_position),by=hgnc_symbol]
matched = match(PathwayGenes$gene, test2$hgnc_symbol)
test2 = test2[matched,]

test3 = all_genes_b38[,unique(end_position),by=hgnc_symbol]
matched = match(PathwayGenes$gene, test3$hgnc_symbol)
test3 = test3[matched,]

test4 = all_genes_b38[,unique(band),by=hgnc_symbol]
matched = match(PathwayGenes$gene, test4$hgnc_symbol)
test4 = test4[matched,]

PathwayGenes[,CHR := test1$V1]
PathwayGenes[,CHR2 := as.numeric(CHR)]
PathwayGenes[is.na(CHR2),CHR2 := 23]
PathwayGenes[,start_position_b38 := test2$V1]
PathwayGenes[,end_position_b38 := test3$V1]
PathwayGenes[,cytoband := paste0(CHR,test4$V1)]
PathwayGenes[CHR %in% c(1:9),cytoband := paste0("0",cytoband)]
setorder(PathwayGenes,CHR2,start_position_b38)

#' # Collapse per cytoband ####
#' ***
PathwayGenes[,table(duplicated(cytoband))]
PathwayGenes[,region_start := start_position_b38 -250000]
PathwayGenes[,region_end := end_position_b38 +250000]

test5 = PathwayGenes[,min(region_start),by=cytoband]
matched = match(PathwayGenes$cytoband,test5$cytoband)
PathwayGenes[,region_start_collapsed := test5[matched,V1]]

test6 = PathwayGenes[,max(region_end),by=cytoband]
matched = match(PathwayGenes$cytoband,test6$cytoband)
PathwayGenes[,region_end_collapsed := test6[matched,V1]]

PathwayGenes[,range := region_end_collapsed - region_start_collapsed]
hist(PathwayGenes$range)

test7 = PathwayGenes[,paste(gene,collapse = " | "),by=cytoband]
matched = match(PathwayGenes$cytoband,test7$cytoband)
PathwayGenes[,genes := test7[matched,V1]]

#' # Save ####
#' ***
save(PathwayGenes, file="../results/03_PathwayGenes_byGene.RData")

PathwayGenes_filt = copy(PathwayGenes)
myNames = c("cytoband", "genes", "CHR", "CHR2", "region_start_collapsed", "region_end_collapsed", "range")
colsOut = setdiff(colnames(PathwayGenes_filt),myNames)
PathwayGenes_filt[,get("colsOut"):=NULL]
setcolorder(PathwayGenes_filt,myNames)
dim(PathwayGenes_filt)
PathwayGenes_filt = PathwayGenes_filt[!duplicated(cytoband),]

#' There is still one overlapping region: 04q13.2 - 04q13.3 
#' 
#' I will merge them manually.
PathwayGenes_filt[cytoband %in% c("04q13.2","04q13.3")]
PathwayGenes_filt[cytoband=="04q13.2",genes := paste0(genes, " | UGT2B4 | UGT2A1 | UGT2A2 | SULT1E1")]
PathwayGenes_filt[cytoband=="04q13.2",region_end_collapsed := 70110145]
PathwayGenes_filt[cytoband=="04q13.2",range := region_end_collapsed - region_start_collapsed]
PathwayGenes_filt[cytoband %in% c("04q13.2","04q13.3")]
PathwayGenes_filt[cytoband=="04q13.2",cytoband := "04q13.2-3"]
PathwayGenes_filt = PathwayGenes_filt[!is.element(cytoband,"04q13.3"),]
PathwayGenes_filt[cytoband %in% c("04q13.2-3")]

save(PathwayGenes_filt, file="../results/03_PathwayGenes_byCytoband.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

