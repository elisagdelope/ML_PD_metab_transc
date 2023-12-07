# Title: lx_shapley.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script calculates shapley values from a given final model and its training set.
# Usage (gene level): Rscript lx_shapley.R 10 (for 10 features)
# Usage (aggregated level): Rscript lx_shapley.R 10 -l pw -s mean 
# Data: data from metabolomics used for training set (only X_train, not Y), ML model.

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
library(fastshap)
library(caret)
library(argparser)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-V0-PD"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_PLOTS <- paste0("../data/", analysis_name , "/01-plots") 
OUT_DIR_PLOTS_SHAP <- paste0(OUT_DIR_PLOTS, "/shap_data") 

# Add command line arguments
target = "DIAGNOSIS"
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "k", help = "number of features selected from DEA", default = 10, short = "k")
p <- add_argument(parser = p, arg = "--level", help = "level of features (metab / pw)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd / pathifier / pca)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
k = as.integer(argv$k) # number of features to select
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level
print(paste("Parameters DEA: k=", k))

if (e_level == "METAB") {
  METAB.FILE <- file.path(OUT_DIR, "log_transformed_V0.tsv")
  IN.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", k, "ft_", target, "_shap_data.csv"))
  IN.MODEL.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_final_model.rds"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", k, "ft_", target, "_shap_values.csv"))
  features_varname <- "METABOLITES_ID"
} else { # aggregations
  METAB.FILE <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, "_V0.tsv"))
  IN.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", k, "ft_", target, "_shap_data.csv"))
  IN.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_final_model.rds"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", k, "ft_", target, "_shap_values.csv"))
  features_varname <- "PATHWAY_NAME"
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(METAB.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
  }
}


# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR_PLOTS_SHAP)) {
    dir.create(OUT_DIR_PLOTS_SHAP, recursive = T)
}



# Data load --------------------------------------------------------------------
shap_data <- vroom(IN.SHAPDATA.FILE, col_types = cols())
ff_model <- readRDS(IN.MODEL.FILE)



# Generate Shapley values ------------------------------------------------------
pfun <- function(object, newdata) {
  predict(ff_model, newdata = newdata, type = "prob")[,2]
}
shap_data <- as.data.frame(shap_data) %>% 
  column_to_rownames(var = "SAMPLE_ID")
shap <- fastshap::explain(ff_model, X = shap_data, pred_wrapper = pfun, nsim = 5)


# output results ---------------------------------------------------------------
readr::write_csv(shap, file = OUT.SHAPVALUES.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()



