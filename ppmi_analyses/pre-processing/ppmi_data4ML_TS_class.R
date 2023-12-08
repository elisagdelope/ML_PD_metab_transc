# Title: ppmi_data4ML_TS_class.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script applies pre-processing steps to generate data for ML modelling of targets PD/UPDRS III from RNAseq longitudinal features (TS), extracted from gene expression and higher-level representations.
# Usage: Rscript ppmi_data4ML_TS_class.R
# Data: data from longitudinal features extracted from transcriptomics data, clinical data (e.g., diagnosis, age, sex).


# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)
library(caret)
library(argparser, quietly = TRUE)
library(matrixStats)




# I/O --------------------------------------------------------------------------
# cmd line arguments
p <- arg_parser("PreprocessCV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--analysis", help = "name of analysis directory (e.g., 02-pred-TS-PD, 02-pred-TS-UPDRS3-class)", required = TRUE)
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
analysis_name <-argv$analysis 
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_DATA <- paste0("../data/", analysis_name , "/05-data4ML")
IN_DIR <- paste0(OUT_DIR, "/temporal_data")
IN_DIR_PATHWAY <- paste0(OUT_DIR_PATHWAY, "/temporal_data")
myseed = 111

if (grepl("TS-PD", analysis_name)) {
  PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
  target = "DIAGNOSIS"
  pheno_4ML <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddlf")) %>%
    rename_with(toupper) %>%
    mutate("PATIENT_ID" = GSEID) %>%
    group_by(PATIENT_ID) %>% 
    dplyr::filter(VISIT == "V08") %>%
    dplyr::select(all_of(c(target, 'PATIENT_ID'))) %>%
    mutate(!!target := case_when(get(target) == "HC" ~ 0,
                                 get(target) == "PD" ~ 1)) 
  
} else if (grepl("TS-UPDRS3", analysis_name)){
  PHENO.FILE <- file.path(OUT_DIR, "ppmi_pheno.tsv")
  target = "UPDRS3_binary"
  pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddlf")) 
  pheno_4ML <- pheno %>%
    dplyr::select(all_of(c(target, 'PATIENT_ID'))) 
}

source("func_data4ML_class.R")



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY)) | (!dir.exists(OUT_DIR_DATA))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
  dir.create(OUT_DIR_DATA, recursive = T)
}

for (temp_feature in c("sd", "lm-time", "lm-lag")) {
  for (e_level in c("GENE", "GOBP", "GOCC", "CORUM")) {
    if (e_level == "GENE") { 
      EXPRS_4ML.FILE <- file.path(IN_DIR, paste(e_level, temp_feature, "temp_expr.tsv", sep = "_"))
      features_varname <- "GENEID"
      process_data4ML_TS(EXPRS_4ML.FILE, pheno_4ML, features_varname, OUT_DIR_DATA, e_level, st, target, temp_feature, myseed, export=TRUE)
    } else { # aggregationss
      
      for (st in c("mean", "median", "sd", "pca", "pathifier")) {
        EXPRS_4ML.FILE <- file.path(IN_DIR_PATHWAY, paste(e_level, st, temp_feature, "temp_expr.tsv", sep = "_"))
        features_varname <- paste0(e_level, "_NAME")
        
        if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(EXPRS.FILE))) { 
          stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
        }
        process_data4ML_TS(EXPRS_4ML.FILE, pheno_4ML, features_varname, OUT_DIR_DATA, e_level, st, target, temp_feature, myseed, export=TRUE)
      }
    }
  }  
}
  


# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


