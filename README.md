# lncGSEA
An R package for linking gene signatures with lncRNA's expression and make prediction of enriched pathways regulated by lncRNAs in human cancer samples

## Prerequisites
```
library(fgsea)
library(data.table)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
```
## Install `lncGSEA` 
```
library(devtools)
install_github("ylab-hi/lncGSEA")
library(lncGSEA)
```
## Download required datasets 
There are two kinds of datasets required for `lncGSEA` to perform its function of finding enriched pathways regulated by lncRNAs.
- lncRNA expression in human cancer samples. 
  - Two public database: mitranscriptome beta and RefLnc. The files were named as "mitranscriptome.expr.fpkm.tsv.gz" (can be downloaded at https://drive.google.com/file/d/15ZucdNxAUT5ZfZZxHBEZ6Q7UvEVjeYfL/view?usp=sharing) and "RefLnc_lncRNA_tumor_sample_FPKM.gz" (download link https://drive.google.com/file/d/1OWyqJlGnN7V0gRh7B-BOIwJSJJ57Zla0/view?usp=sharing), respectively. 
- gene expression matrix for each cohort in TCGA study.
  - Example: PRAD.FPKM.txt, BRCA.FPKM.txt, COAD.FPKM.txt

All datasets can be downloaded from this shared link:

## Create a data folder in your current working directory
Please create a data folder by the following command to store the downloaded datasets.

```
if (!file.exists("data")){
    dir.create("data")
}
```

## Examples

Create an expression data frame for PRAD cohort, columns are one of transcripts of ARLNC1 (e.g. ENST00000561519) and other genes. Rows are tumor samples from PRAD cohort. The first column name must be the cohort name one use in `pre_gsea` and the second column name is the transcript id. A same transcript id may have different versions, `pre_gsea` use transcript id without version number ("ENST00000561519" instead of "ENST00000561519.[0-9]"). 

```
test <- pre_gsea("PRAD", "ENST00000561519")
test[1:4, 1:4]
```

The first 4 rows and columns are shown below:
    


                                  PRAD ENST00000561519.5 TSPAN6  DPM1 SCYL3
     1754 TCGA-G9-7521-01A-11R-2263-07            0.0104  10.62 24.92 1.971
     1757 TCGA-KK-A7AU-01A-11R-A32O-07            0.2343  13.35 23.39 3.539
     1760 TCGA-EJ-7125-01A-11R-1965-07            0.7116  18.71 15.37 1.835
     1763 TCGA-ZG-A9LM-01A-11R-A41O-07            0.4569  16.89 29.95 2.484
     1765 TCGA-QU-A6IM-01A-11R-A31N-07            0.0000   8.45  9.04 0.789
