# Title: ppmi_generate_temp_features.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates temporal features from gene expression data
# Usage (gene level): Rscript ppmi_generate_temp_features.R 
# Usage (aggregated level): Rscript ppmi_generate_temp_features.R -l gobp -s mean
# Data: data from expression counts at several timepoints (>=3), clinical data (diagnosis).

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)
library(argparser, quietly = TRUE)
library(matrixStats)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-TS-PD"
IN_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
IN_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR <- paste0(IN_DIR, "/temporal_data")
OUT_DIR_PATHWAY <- paste0(IN_DIR_PATHWAY, "/temporal_data")
#IN_DIR <- paste0("../", analysis_name , "/02-outfiles") 
#IN_DIR_PATHWAY <- paste0("../", analysis_name , "/04-pathway_level")

# Add command line arguments
p <- arg_parser("Generate temporal features from gene expression data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level

if (e_level == "GENE") { 
  EXPRS.FILE <- file.path(IN_DIR, "flt_norm_star_all_TS.tsv")
  features_varname <- "GENEID"
} else {  # aggregations
  EXPRS.FILE <- file.path(IN_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
  features_varname <- paste0(e_level, "_NAME")
  if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(EXPRS.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_generate_temp_features.R --help for the right usage.")
  }
}



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------------
# star raw counts & (pre-filtered) clinical data.
expr <- vroom(EXPRS.FILE, col_types = cols()) %>% 
  rename_with(toupper)
expr = expr[!duplicated(expr[[features_varname]]),]



# Data reformatting ------------------------------------------------------------
# remove patients that DON'T have data for 4 timepoints
n_tp = 4
patients_list <- expr %>%
  pivot_longer(cols = stringr::str_subset(colnames(expr), "[0-9]{4}\\.[A-Z]"), names_to = "SAMPLE", values_to = "EXPRS") %>%
  separate(SAMPLE, c("PATIENT_ID", "VISIT"), sep = '\\.') %>%
  group_by(!!as.symbol(features_varname), PATIENT_ID) %>% 
  filter(n() == n_tp) %>%   
  ungroup() %>% 
  dplyr::select(PATIENT_ID) %>%
  distinct() %>%
  pull()

expr <- expr %>% 
  select(matches(c(features_varname, patients_list)))

# calculate dynamic features from expression
source("func_temporal_features.R")
for (temp_feature in c("sd", "lm-time", "lm-lag")) {
  if (temp_feature == "sd") {
    expr_4ML <- exprs_sd(expr, features_varname)
  } else if (temp_feature == "lm-time") {
    expr_4ML <- exprs_lm(expr, features_varname)
  } else if (temp_feature == "lm-lag") {
    expr_4ML <- exprs_lm_lag(expr, features_varname)
  } 
  
  # remove those genes having NA for the dynamic features -----
  expr_4ML <- expr_4ML %>%
    select_if(~ !any(is.na(.))) 
  
  # output results ------
  if (e_level == "GENE") { 
    OUT.FILE <- file.path(OUT_DIR, paste(e_level, temp_feature, "temp_expr.tsv", sep = "_"))
  } else { 
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, temp_feature, "temp_expr.tsv", sep = "_"))
  }
  readr::write_tsv(expr_4ML, file = OUT.FILE)
  
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()



