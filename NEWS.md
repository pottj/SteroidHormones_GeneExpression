# NEWS file

Here I want to track my progress

## August, 21st, 2025

- update TWAS with replication in SHIP
- update sample selection & checks in UKB and rerun PWAS 

## May, 9th, 2025

- rearrange and update repository: 
  - one folder for LIFE analyses
  - one folder for UKB analyses
  - one folder for SHIP replication

## December, 13th, 2024

- MR-Egger: 3 genes (*EIF4A1P2*, *EIF5AL1*, *CNST*)
- MR-IVW with pathway SNPs: 1 gene (*EIF4A1P2*)
- MR-Egger with pathway SNPs: 2 genes (*EIF4A1P2*, *EIF5AL1*)
- MVMR-IVW with TT pathway SNPs & all BMI SNPs: 6 genes (*MYLIP*, *ABCA1*, *EIF5AL1*, *EIF4A1P2*, *CITED4*, *ABCG1*)   

## Novmeber, 29th, 2024

- 2SLS: no significant results after hier. FDR (transcriptome-wide and gene set based)
- 1-sample MR using PGS as instruments: after hier. FDR, *SERPINA1* is still significantly linked to CORT. But this might just be confounding by LD (CORT PGS uses only SNPs in the *SERPINA6* gene region, neighbor of *SERPINA1*)
- 2-sample MR using SNPs as instruments: 
    - data source for trans-eQTLS: eQTLGen phase I 
    - only TT and E2, because after matching no significant SNP remained for CORT
    - after hier. FDR, three genes are significantly linked to TT in men :
        - *ABCA1*
	      - *ABCG1*
	      - *EIF4A1P2*
- To do:  
    - rerun 2-sample MR using MR-Egger (high pleiotropy)
    - rerun 2-sample MR using only pathway SNPs (at the moment all ~100 SNPs that are significantly associated with TT and in eQTLGen)
- Discussed with Steve: 
    - use MVMR to adjust for BMI effect
    - maybe worthwhile to select genes with known AR binding site, and with either known BMI association or no BMI link whatsoever. Then analyze the BMI - TT feedback loop in more detail
    - is it BMI or is it fat mass / adipose tissue? (but we do not have adipose tissue gene expression, GTEx only provides trans-eQTLs on the chromosome of the respective gene)
    
## November, 22nd, 2024

**PGS**:

- update: choose SNPs from the same set: first select SNPs for Adult (pruning, clumping), then use these SNPs for Heart as well (better for comparison)
- PGS selection: (based on criteria genome-wide significant and highest $r^2$ per hormone and sample set in LIFE-Adult)

| hormone | sample set | p-value threshold | SNP set      | comment                     |
| ------- | ---------- | ----------------: | ------------ | --------------------------- |
| TT      | men        | 0.001             | genome-wide  | pathway-wide as sensitivity | 
| TT      | women      | 0.001             | genome-wide  | pathway-wide as sensitivity | 
| E2      | men        | 0.001             | pathway-wide | genome-wide as sensitivity  | 
| E2      | women      |                   |              | no strong instruments       | 
| CORT    | combined   | 0.001             | *SERPINA6*   | p-value<1e-4 as sensitivity | 

**TWAS**

- rerun restricted to relevant traits
    - TT (separated by sex), E2 (separated by sex), and CORT (sex-combined)
    - hierarchical FDR over the 5 traits and all ~20,000 probes (a bit conservative, as women each probe-individual pair is only tested 3 times, but it is the easiest way)
    - two models: with and without BMI adjustment - check hierarchical FDR over 10 trait - model combinations
- check overlap between 
    - the two models in each study (venn2 + euler plot to check effect direction)
    - LIFE Adult and LIFE Heart (venn4 + euler plot to check effect direction)
- create bar plots and volcano plots

Result in LIFE Adult: (hierarchical FDR significant probes) 

| hormone | sample set | # probes (no BMI) | # probes (BMI adj.) | # probes in overlap | up  | down | 
| ------- | ---------- | ----------------: | ------------------: | ------------------: | ---:| ---: |
| TT      | men        | 7657              | 281                 | 280                 | 110 | 170  |
| E2      | men        | 10                | 10                  | 5                   | 0   | 5    |
| TT      | women      | 3                 | 2                   | 0                   | 1   | 1    |
| E2      | women      | 1473              | 77                  | 77                  | 18  | 59   |
| CORT    | combined   | 973               | 4421                | 930                 | 415 | 515  |

Replicated in LIFE Heart: (p-value < 0.05 and same effect direction in Heart)

| hormone | sample set | # probes replicated | up  | down | different direction | 
| ------- | ---------- | ------------------: | --: | ---: | -------------------:| 
| TT      | men        | 53                  | 10  | 40   | 3                   |
| E2      | men        | 0                   |     |      |                     |
| TT      | women      | 0                   |     |      |                     |
| E2      | women      | 2                   | 1   |      | 1                   |
| CORT    | combined   | 565                 | 276 | 275  | 14                  |

## November, 15th, 2024

**PGS**:

- updated the scores using genome-wide SNPs, we only have **strong instruments for TT** in men and women 
- rerun scores using SNPs within/near steroid hormone pathway enzymes encoding genes
- run PGS special for CORT: different p-value thresholds, using gene dosage matrix not in PLINK but in R
- To do: final selection of score per hormone, using LIFE Adult as main analysis

## August 2024

- 02/08/2024: I decided to rerun a candidate GWAS for all hormones to use the same genetic data set throughout the project
    - update data prep again
    - extract SNPs and run GWAS
    - PROBLEM with LIFE-Heart - mail to Katrin & Holger 
    - some tests for the TWAS and 2SLS

## July 2024

- 05/07/2024: 
    - updated data prep for LIFE Heart (data download successful on 02/07/2024)
    - updated data prep for LIFE Adult
    - TWAS and hierarchical FDR in Adult 
    - TWAS and hierarchical FDR in Heart 
    - TWAS ToDo: get genes associated in both Adult and Heart, get Vulcano Plots
    - downloaded GWAS summary statistics and prepped data for PGS

### Data prep

So far so good. Might want to discuss the decisions for premenopausal women, but the rest should be okay. Maybe I will just evaluate the men and post-menopausal women. 

Check E2: There were about 200 post-menopausal women in LIFE-Adult, who had their ovaries/uterus removed. Given that they are postmenopausal, this might not be critical, and removal can be used as a covariate in the TWAS model. But I did not want to lose all those women. 

### PGS 

Not yet completely sure how the best procedure is. First step is to get the regions of genes coding for steroid hormone pathway genes (using KEGG and biomaRt ensembl). Second step is to extract SNPs in that regions (publicly available summary statistics). Third is to create a score using the effect sizes as weights (different combinations of SNPs / genes). Finally, I check if the scores are associated with the hormones and if they are good instruments. If not, then I have to run a candidate GWAS, and select the best SNPs in LIFE alone... 

Data sources: 

- UKB TT (Ruth et al.), hg38 
- E2 (Schmitz et al.), binary trait (<LOD), hg38 
- CORT (Chan and Wu), CORNET + LIFE, hg37 liftover using [Bioconductor](https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP150.GRCh38.html)
- LIFE-Adult and LIFE-Heart TopMed imputation (hg38)

### TWAS

Fixed matching problem by using the sample annotation table from the QC data.

## June 2024

- 14/06/2024: 
    - create github repository
    - create first script for data preparation of LIFE-Adult
- 21/06/2024: not worked on this project (participation at the MR Conference in Bristol)
- 28/06/2024:
    - mail with Yvonne Dietz (LIFE-DM): 
        - steroid hormone data for LIFE-Heart is missing / not in the shared nextCloud repository
        - Yvonne added the data, but I cannot download them
        - possible explanation: admins update ISF2 
        - try again on Monday
    - finished data prep for LIFE-Adult (after applying all sample exclusion criteria): 
        - **TWAS**: There are 2183 samples with hormone AND GE data
        - **PGS**:  There are 4535 samples with hormone AND genetic data 
        - **TSLS**: There are 2064 samples with hormone AND genetic AND GE data
        - **eQTL**: There are 2301 samples with genetic AND GE data 
    - finished data prep for LIFE-Heart (after applying all sample exclusion criteria; **this is done using data from a previous PV (SASHA, 505), and I will update the data as soon as possible, but it should not change much**):
        - **TWAS**: There are 1477 samples with hormone AND GE data
        - **PGS**:  There are 1602 samples with hormone AND genetic data 
        - **TSLS**: There are 1477 samples with hormone AND genetic AND GE data
        - **eQTL**: There are 2483 samples with genetic AND GE data 
    - started **TWAS** for LIFE-Adult
        - stopped because matching problem between D00403 and GE data - continue here!
  
## April/May 2024

Write and submit proposal for LIFE Project Agreements (LIFE PV)

- 23/04/2024: proposal sent out to collaborators for approval
- 29/04/2024: proposal submitted to LIFE committee
- 13/05/2024: LIFE-PV meeting to discuss proposal; approval on condition that all individual-level data remain in Leipzig, and all computations will be performed on the IMISE servers
- 27/05/2024: data release; **PV number 785**