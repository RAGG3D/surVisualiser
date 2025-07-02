install.packages("BiocManager")
BiocManager::install("tidybulk")
BiocManager::install("tidyverse")
BiocManager::install("tidyHeatmap")
# BiocManager::install("clusterProfiler")
# BiocManager::install("illuminaHumanv4.db")
# BiocManager::install("preprocessCore")
# BiocManager::install("singscore")
# BiocManager::install("singscore")
# BiocManager::install("textshaping")
# BiocManager::install("WGCNA")
options(connectionObserver = NULL)
devtools::install_github("stemangiola/tidybulk")
library(cowplot)
library(stringr)
library(vegan)
library(ggplotify)
library(PCAtools)
library(forcats)
library(TCGAbiolinks)
library(WGCNA)
library(Polychrome)
library(igraph)
library(ggraph)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(org.Hs.eg.db)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidybulk)
library(tidyHeatmap)
library(tidymodels)
library(data.table)
library(survminer)
library(survival)
library(cowplot)
library(scales)
library(ggsci)
library(foreach)
library(reshape2)
library(Hmisc)
library(gridExtra)
library(formattable)
library(psych)
library(clusterProfiler)
library(ggpmisc)
library(ggrepel)
library(Seurat)
library(patchwork)
library(illuminaHumanv4.db)
library(preprocessCore)
library(grid)
library(ggsci)
library("RColorBrewer")
library(ggeasy)
library(singscore)
library(GSEABase)
library(hacksig)
library(random)
library(msigdbr)
library(ggforce)

# #url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"
# #download.file(url, destfile = "kegg_2021.gmt")
# kegg <- read.gmt("kegg_2021.gmt")
# 
# #url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human"
# #download.file(url, destfile = "wp.gmt")
# wp <- read.gmt("wp.gmt")
# 
# #url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021"
# #download.file(url, destfile = "gobp.gmt")
# gobp <- read.gmt("gobp.gmt")
# 
# #url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019"
# #download.file(url, destfile = "biop.gmt")
# biop <- read.gmt("biop.gmt")
# 
# #url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016"
# #download.file(url, destfile = "react.gmt")
# react <- read.gmt("react.gmt")


NKT = c('ReNK', 'IL2NK', 'SPANK', "T Helper", "Naive CD8 T", "GD T", 'CD4 Tcm', 
        'CD4 Tem', 'CD8 Tcm',  'CD8 Tem', "Treg")

#*Aggregate samples and do TMM normalization
TCGA_transcript <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}

TCGA_transcript_tumor <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/tumor_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}


TCGA_transcript_normal <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/normal_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}


TCGA_transcript_tumor_raw <- function(cancer){
  map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
          function(i){
            data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
              mutate(sample = i)
          }) %>%
    tidybulk::rename(raw_count = `V2`) %>%
    mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
    inner_join(toTable(org.Hs.egENSEMBL)) %>%
    inner_join(toTable(org.Hs.egSYMBOL)) %>%
    dplyr::select(sample, symbol, raw_count) %>%
    mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
    inner_join(
      read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/tumor_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
    mutate(sample = `Case ID`) %>%
    dplyr::select(sample, symbol, raw_count)
  
}

#Combine clinical data
clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("/home/yzsun/TCGA clinical data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}

