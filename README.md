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
