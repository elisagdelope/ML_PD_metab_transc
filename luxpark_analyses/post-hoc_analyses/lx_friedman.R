# Title: lx_friedman.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs comparisons among model performance, data types performance, and feature selection methods.
# Usage: Rscript lx_friedman.R 
# Usage (can be run from command line): Rscript lx_friedman.R > ../reports/PRED-V0-TS-PD_results/stdout_friedman.txt 2>err/stderr_friedmans.txt
# Data: data from summary of results xlsx file.


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
library(PMCMRplus)
library(readxl)



# I/O --------------------------------------------------------------------------
IN_DIR <- "../reports/PRED-V0-TS-PD_results"
IN_FILE <- file.path(IN_DIR, "summary_AUC_scores.tsv")
CONFOUNDERS_FILE <- file.path(IN_DIR, "CONFOUNDERS_summary_results.csv")
ORIG_FILE <- file.path(IN_DIR, "METAB_summary_results.xlsx")
OUT_DIR <- IN_DIR
OUT_DIR_PLOTS <- "../reports/PRED-V0-TS-PD_results/plots"
source("func_friedman_connover_holm_nxn.R")
source("func_friedman.R")


# Data transformations ---------------------------------------------------------

if (file.exists(IN_FILE)){
  auc <- vroom(IN_FILE, show_col_types = FALSE)
} else {
  # Data transformations ---------------------------------------------------------
  data <- list()
  data[["lasso_ftsel"]] <- read_excel(ORIG_FILE, sheet = 1)
  data[["limma_ftsel"]] <- read_excel(ORIG_FILE, sheet = 2)
  data[["rf_noftsel"]] <- read_excel(ORIG_FILE, sheet = 3)
  
  data_types <- unique(data[["lasso_ftsel"]][,1])
  data_types <- data_types[(!is.na(data_types)) & data_types != "de novo"]
  
  data[["lasso_ftsel"]] <- data[["lasso_ftsel"]] %>%
    filter(.[[1]] != "de novo") %>%
    rename_at(1,~"data_type") %>%
    dplyr::select(-2) %>%
    mutate(across(-1, as.double)) %>% 
    mutate(linearSVM=as.double(linearSVM))
  
  data[["limma_ftsel_10"]] <- data[["limma_ftsel"]] %>%
    filter(`N features` == 10) %>%
    rename_at(1,~"data_type") %>%
    mutate("data_type" = data_types) %>%
    dplyr::select(-c(2,3)) %>%
    mutate(across(-1, as.double)) %>% 
    dplyr::select(-"LogitBoost")
  
  data[["limma_ftsel"]] <- data[["limma_ftsel"]] %>%
    filter(`N features` == 100) %>%
    rename_at(1,~"data_type") %>%
    mutate("data_type" = data_types) %>%
    dplyr::select(-c(2,3)) %>%
    mutate(across(-1, as.double)) %>% 
    dplyr::select(-"LogitBoost")
  
  data[["rf_noftsel"]] <- data[["rf_noftsel"]] %>% 
    filter(!is.na(.[[1]])) %>%
    rename_at(1,~"data_type") %>%
    dplyr::select(-2) %>%
    mutate(across(-1, as.double)) 
  
  an_datatype <- list()
  an_model <- list()
  type_model_auc <- list()
  for (i in names(data)) {
    
    an_datatype[[i]] <- as.data.frame(apply(data[[i]] %>% 
                                              column_to_rownames(var = "data_type"), 1, median))
    colnames(an_datatype[[i]]) <- "AUC"
    an_datatype[[i]] <- an_datatype[[i]] %>% 
      rownames_to_column(var = "data_type")
    
    an_model[[i]] <- as.data.frame(apply(data[[i]] %>% 
                                           column_to_rownames(var = "data_type") %>% t, 1, median))
    colnames(an_model[[i]]) <- "AUC"
    an_model[[i]] <- an_model[[i]] %>% 
      rownames_to_column(var = "model")
    
    type_model_auc[[i]] <- data[[i]] %>% 
      pivot_longer(cols = -data_type, names_to = "model", values_to = "AUC")
    
  }
  
  
  
  # export AUC scores table ------------------------------------------------------
  #  merge all AUC scores into full AUC table
  auc <- bind_rows(type_model_auc, .id = "ftsel")
  # print full AUC table
  readr::write_tsv(auc, file = file.path(OUT_DIR, "summary_AUC_scores.tsv"))
}


# RF Lasso vs. other feature selection -----------------------------------------
print("################## FT. SELECTION COMPARISON ANALYSIS ##################")

type_model_auc_RF <- auc %>%
  dplyr::filter(model=="RandomForest" | model=="rf")  %>%
  mutate(model = case_when(model == "rf" ~ "RandomForest", TRUE ~ model))

print(paste("methods compared:", paste(unique(type_model_auc_RF$ftsel), collapse=", ")))

frd <- friedman(groups=type_model_auc_RF$ftsel, 
                blocks=type_model_auc_RF$data_type,
                score=type_model_auc_RF$AUC)
if (frd$p.value < 0.05) {
  fr_nx1 <- friedman_holm_nx1(groups=type_model_auc_RF$ftsel, 
                              blocks=type_model_auc_RF$data_type,
                              score=type_model_auc_RF$AUC,
                              col_control=1)
  print(paste("Friedman + Holm post-hoc test: ", names(which.min(fr_nx1$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(fr_nx1$test[2,][which(fr_nx1$test[2,] < 0.05)]), collapse=", ")))
} else{
  print("Friedman hypothesis test was not rejected")
}


# Pairwise comparison: models x models -----------------------------------------
print("################## MODELS COMPARISON ANALYSIS ##################")
type_model_auc_i <- auc %>%
  dplyr::filter(ftsel=="lasso_ftsel")
frd <- friedman(groups=type_model_auc_i$model, 
                blocks=type_model_auc_i$data_type,
                score=type_model_auc_i$AUC)
if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_i$model,  
                                        blocks=type_model_auc_i$data_type,
                                        score=type_model_auc_i$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_models.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
}


# Pairwise comparison: data types x data types ---------------------------------
print("################## FEATURE TYPES COMPARISON ANALYSIS ##################")
type_model_auc_i <- auc %>%
  dplyr::filter(ftsel=="lasso_ftsel")
frd <- friedman(groups=type_model_auc_i$data_type,
                blocks=type_model_auc_i$model,
                score=type_model_auc_i$AUC)
if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_i$data_type,
                                        blocks=type_model_auc_i$model,
                                        score=type_model_auc_i$AUC)
  
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_datatypes.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
  
} else{
  print("Friedman hypothesis test was not rejected")
}



# Pairwise comparison: CONFOUNDERS x CONFOUNDERS ---------------------------------
print("################## CONFOUNDERS COMPARISON ANALYSIS ##################")

conf_auc <- vroom(CONFOUNDERS_FILE, show_col_types = FALSE)
colnames(conf_auc)[1] <- "vars_type"
conf_auc <- conf_auc %>% 
  pivot_longer(cols = -vars_type, names_to = "model", values_to = "AUC")

frd <- friedman(groups=conf_auc$vars_type,
                blocks=conf_auc$model,
                score=conf_auc$AUC)
if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=conf_auc$vars_type,
                                        blocks=conf_auc$model,
                                        score=conf_auc$AUC)
  
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  
  # print plot
  fr_nxn$plot <- fr_nxn$plot  + labs(x = "Types of covariates", y = "Types of covariates") + theme(
    axis.text.x = element_text(size = 12, angle = 14, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),    # Set x-axis label size
    axis.title.y = element_text(size = 14),     # Set y-axis label size
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_confounders.pdf"), width = 9, height = 4)
  print(fr_nxn$plot)
  dev.off()
  
} else{
  print("Friedman hypothesis test was not rejected")
}

fr_nx1 <- friedman_holm_nx1(groups=conf_auc$vars_type, 
                            blocks=conf_auc$model,
                            score=conf_auc$AUC,
                            col_control=5)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

