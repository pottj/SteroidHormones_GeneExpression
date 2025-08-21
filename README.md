# Mendelian Randomization of steroid hormones on gene expression in whole blood

created: 14/06/2024

last updated: 10/05/2024

## Overview 

This project is a collaboration of me, Stephen Burgess (MRC BSU, University of Cambridge, UK), Markus Scholz (IMISE, University of Leipzig, GER), and Alexander Teumer (University of Greifswald, GER) to test if we can detect causal effects of steroid hormones on gene expression levels in whole blood / PBMCs. 

At the moment, we are focusing on gene expression using individual-level data from the LIFE-Adult and LIFE-Heart studies (Leipzig, GER; project agreement 785), and from SHIP (Greifswald, GER; project agreement SHIP/2025/31/D), and proteomics data from the UK Biobank (UKB, UK; project application 98032). 

## Short analysis plan

### Gene expression analyses (LIFE and SHIP)

1) Data preparation

Check and filter data, create groups (post-/peri-menopausal), check hormonal distribution and decide for or against transformation. 

2) Polygenetic Score (PGS)

Generation of scores using publicly available data (PGS Catalog, previously published GWAS summary statistics); testing the scores for validity in MR analyses by linear regression models. 

3) Transcriptome-wide association studies (TWASs)

Testing the association of steroid hormones on gene expression using limma. Results will be corrected for multiple testing using hierarchical FDR (correcting for multiple tests on genes and hormones). For this step, we will use all individuals of LIFE-Adult and LIFE-Heart with gene expression data and available steroid hormone measurements. Genes associated in both LIFE-Adult and LIFE-Heart will be selected for further analyses (candidate genes).

4) Mendelian Randomization (MR)

Testing with 2-stage-least-square for causal effect of hormones on gene expression in whole blood (1-sample approach). In addition, testing with IVW and Egger for causal effects using data from the UKB and eQTLGen. 

5) Multivariable MR (MVMR)

Testing for independent hormonal effects, conditioned on BMI effect. 

### Proteomics analyses (UK Biobank)

1) Data preparation 
2) Proteome-wide association study (PWAS)
3) MR using 2SLS 

### Data used so far

**GWAS summary statistics**

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
  
**Used gene regions**: 

- 62 genes as listed in the [KEGG Steroid hormone biosynthesis pathway](https://www.genome.jp/dbget-bin/www_bget?path:hsa00140)
- *SHBG* (Sex Hormone Binding Globulin)
- *SERPINA6* (aka CBG, Corticosteroid Binding Globulin)

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
