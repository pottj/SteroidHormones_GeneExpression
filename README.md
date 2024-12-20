# Mendelian Randomization of steroid hormones on gene expression in whole blood

created: 14/06/2024

last updated: 04/07/2024

## Overview 

This project is a collaboration of me, Stephen Burgess (MRC BSU, University of Cambridge, UK), and Markus Scholz (IMISE, University of Leipzig, GER) to test if we can detect causal effects of steroid hormones on gene expression levels in whole blood / PBMCs. At the moment, we are using individual-level data from the LIFE-Adult and LIFE-Heart studies (Leipzig, GER; project agreement 785). 

## Short analysis plan

### PGS

Generation of scores using publicly available data (PGS Catalog, previously published GWAS summary statistics); testing the scores for validity in MR analyses by linear regression models. For this step, we will use all individuals of LIFE-Adult and LIFE-Heart with genetic data and available steroid hormone measurements.

Used data sets (so far): 

- Testosterone (TT): [Ruth et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32042192/)
    - [GWAS Catalog GCST90012112](https://www.ebi.ac.uk/gwas/studies/GCST90012112) (women)
    - [GWAS Catalog GCST90012113](https://www.ebi.ac.uk/gwas/studies/GCST90012113) (men)
- Estradiol (E2): [Schmitz et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34255042/), binary trait
    - [GWAS Catalog GCST90020092](https://www.ebi.ac.uk/gwas/studies/GCST90020092) (women)
    - [GWAS Catalog GCST90020091](https://www.ebi.ac.uk/gwas/studies/GCST90020091) (men)
- Cortisol (CORT): [Chan and Wu](https://pubmed.ncbi.nlm.nih.gov/38525495/), (CORNET consortium + LIFE)
    - [figshare](https://figshare.com/articles/dataset/cortisol_cornet_life_combined/26182004) 
- eQTLGen (trans-eQTLs): [Võsa, U., Claringbould, A., (…), Franke, L.; 2018; Unraveling the polygenic architecture of complex traits using blood eQTL meta-analysis](https://eqtlgen.org/publications.html), [publication](https://pubmed.ncbi.nlm.nih.gov/34475573/)
    - [Full trans-eQTL summary statistics download](https://eqtlgen.org/trans-eqtls.html)
  
  
Used gene regions: 

- 62 genes as listed in the [KEGG Steroid hormone biosynthesis pathway](https://www.genome.jp/dbget-bin/www_bget?path:hsa00140)
- *SHBG* (Sex Hormone Binding Globulin)
- *SERPINA6* (aka CBG, Corticosteroid Binding Globulin)

Used hormones: 

- Testosterone (TT, sex-stratified)
- Estradiol (E2, sex-stratified)
- Progesterone (P4, sex-stratified)
- Cortisol (CORT, sex-combined)
- Aldosterone (ALDO, sex-combined)

Used SNPs: 

- PGS SNPs
- suggestive significant SNPs (pairwise independent)
- SNPs within/near genes coding for enzyms of the steroid hormone biosynthesis pathway

### TWAS

Testing the association of steroid hormones on gene expression using limma. Results will be corrected for multiple testing using hierarchical FDR (correcting for multiple tests on genes and hormones). For this step, we will use all individuals of LIFE-Adult and LIFE-Heart with gene expression data and available steroid hormone measurements. Genes associated in both LIFE-Adult and LIFE-Heart will be selected for further analyses (candidate genes).

### 2SLS / 2SMR / SEM

Testing for causal effects of steroid hormones on gene expression of selected candidate genes. For this step, we will use all individuals of LIFE-Adult and LIFE-Heart with genetic data (PGS), gene expression data, and available steroid hormone measurements. Genes significantly affected in both LIFE-Adult and LIFE-Heart and with the same effect direction in both studies will be selected for further enrichment analyses.

### Enrichment

We will test the selected genes for pathway enrichment (e.g. pathways obtained from KEGG), enrichment of hormone response elements (consensus sequence), and for overlaps between the strata (shared causal effects between men and women or shared causal effects between pre- and postmenopausal women). These analyses will be done using the summary statistics from the TSLS regression.

## List of abbreviations

| Abbreviation | Meaning                                 | 
| :----------- | :-------------------------------------- | 
| SH           | Steroid hormone                         |
| GE           | Gene expresssion                        |
| TT           | Total testosterone                      |
| E2           | Estradiol                               |
| P4           | Progesterone                            |
| CORT         | Cortisol                                |
| ALDO         | Aldosterone                             |
| PGS          | Polygenetic Score                       |
| TWAS         | Transcriptome-wide association analysis |
| MR           | Mendelian Randomization                 |
| SNP          | Single nucleotide polymorphism          |
| FDR          | False discovery rate                    |
| 2SLS         | 2-stage least square regression         |
| 2SMR         | 2-sample Mendelian Randomization        |
| SEM          | Structure equation model                |
