# NEWS file

Here I want to track my progress

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