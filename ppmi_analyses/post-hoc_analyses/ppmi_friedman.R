# Title: ppmi_friedman.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs comparisons among model performance, data types performance, and feature selection methods built for PD/control sample classification.
# Usage: Rscript ppmi_friedman.R 
# Usage (can be run from command line): Rscript ppmi_friedman.R > ../reports/PRED-BL-TS-PD_results/stdout_friedman.txt 2>err/stderr_friedman.txt
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
source("func_friedman_connover_holm_nxn.R")
source("func_friedman.R")

IN_DIR <- "../reports/PRED-BL-TS-PD_results"
IN_FILE_BL <- file.path(IN_DIR, "summary_AUC_scores_BL.tsv")
IN_FILE_BLTS <- file.path(IN_DIR, "summary_AUC_scores_BLTS.tsv")
IN_DIR_BL = "../02-pred-BL-PD/02-outfiles"
IN_DIR_PATHWAY_BL = "../02-pred-BL-PD/04-pathway_level"
IN_DIR_TS = "../02-pred-TS-PD/02-outfiles"
IN_DIR_PATHWAY_TS = "../02-pred-TS-PD/04-pathway_level"
OUT_DIR <- IN_DIR
OUT_DIR_PLOTS <- "../reports/PRED-BL-TS-PD_results/plots"
OUT_DIR_BL <- IN_DIR_BL
OUT_DIR_TS <- IN_DIR_TS
PW.FILES <- list()
GENE.FILES <- list()
target = "DIAGNOSIS"
ORIG_FILE_BL_DESEQ <- file.path(IN_DIR, "AUC_analysis_BL.xlsx")
ORIG_FILE_LASSO <- file.path(IN_DIR, "Tables AUC-analyses_manual.xlsx")


# Analyses at BL ---------------------------------------------------------------

if (file.exists(IN_FILE_BL)){
  auc <- vroom(IN_FILE_BL, show_col_types = FALSE)
} else {
  # Data transformations -------------------------------------------------------
  data <- list()
  data[["deseq"]] <- read_excel(ORIG_FILE_BL_DESEQ, sheet = 1)
  data[["rf_corr"]] <- read_excel(ORIG_FILE_BL_DESEQ, sheet = 2)
  data[["rf_noftsel"]] <- read_excel(ORIG_FILE_BL_DESEQ, sheet = 3)
  data[["lasso"]] <- read_excel(ORIG_FILE_LASSO, sheet = 8)
#  data_types <- unique(data[["deseq"]][,1])
#  data_types <- data_types[(!is.na(data_types)) & data_types != "de novo"]
  
  data[["deseq_10"]] <- data[["deseq"]] %>%
    rename_at(1,~"data_type") %>%
    dplyr::select(-contains("Final"), -contains("RF_")) %>%
    mutate(across(-1, as.double)) %>%
    dplyr::filter(grepl("_10", data_type) & !grepl("_100|_1000", data_type)) %>%
    mutate(data_type = gsub("_10", "", data_type))
  
  data[["deseq_100"]] <- data[["deseq"]] %>%
    rename_at(1,~"data_type") %>%
    dplyr::select(-contains("Final"), -contains("RF_")) %>%
    mutate(across(-1, as.double)) %>%
    dplyr::filter(grepl("_100", data_type) & !grepl("_1000", data_type)) %>%
    mutate(data_type = gsub("_100", "", data_type))
  
  # deseq 1000 filter will not be used as there were <1000 CORUM protein complexes -> no ft selection for those.
  #data[["deseq_1000"]] <- data[["deseq"]] %>%
  #  rename_at(1,~"data_type") %>%
  #  dplyr::select(-contains("Final"), -contains("RF_")) %>%
  #  mutate(across(-1, as.double)) %>%
  #  dplyr::filter(grepl("_1000", data_type)) %>%
  #  mutate(data_type = gsub("_1000", "", data_type))
  
  data[["deseq"]] <- NULL

  data[["rf_corr"]] <- data[["rf_corr"]] %>% 
    rename_at(1,~"data_type") %>%
    unite(data_type, data_type, st, sep = "_", remove = TRUE)  %>%
    mutate(data_type = gsub("GENE_NA", "GENE", data_type)) %>%
    mutate(across(-1, as.double)) 
  
  data[["rf_noftsel"]] <- data[["rf_noftsel"]] %>% 
    rename_at(1,~"data_type") %>%
    unite(data_type, data_type, st, sep = "_", remove = TRUE)  %>%
    mutate(data_type = gsub("GENE_NA", "GENE", data_type)) %>%
    mutate(across(-1, as.double)) 
  
  data[["lasso"]] <- data[["lasso"]] %>%
    fill(1) %>%
    rename_at(1,~"data_type") %>%
    drop_na(2) %>%
    unite(data_type, data_type, 2, sep = "_", remove = TRUE)  %>%
    mutate(data_type = gsub("GENE_-", "GENE", data_type)) %>%
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
  readr::write_tsv(auc, file = file.path(OUT_DIR, "summary_AUC_scores_BL.tsv"))
}



# RF Lasso vs. other feature selection at BL -----------------------------------
print("################## FT. SELECTION COMPARISON ANALYSIS (BL) ##################")

type_model_auc_RF <- auc %>%
  dplyr::filter(model=="RandomForest" | model=="rf")  %>%
  mutate(model = case_when(model == "rf" ~ "RandomForest", TRUE ~ model)) %>%
  dplyr::filter(grepl("GENE|mean|sd", data_type)) # only mean, sd and GENE will be used for comparison of Lasso.

print(paste("methods compared:", paste(unique(type_model_auc_RF$ftsel), collapse=", ")))

frd <- friedman(groups=type_model_auc_RF$ftsel, 
                blocks=type_model_auc_RF$data_type,
                score=type_model_auc_RF$AUC)
if (frd$p.value < 0.05) {
  fr_nx1 <- friedman_holm_nx1(groups=type_model_auc_RF$ftsel, 
                              blocks=type_model_auc_RF$data_type,
                              score=type_model_auc_RF$AUC,
                              col_control=3)
  print(paste("Friedman + Holm post-hoc test: ", names(which.min(fr_nx1$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(fr_nx1$test[2,][which(fr_nx1$test[2,] < 0.05)]), collapse=", ")))
} else{
  print("Friedman hypothesis test was not rejected")
}

# Lasso vs. DEGs with RF, Log. R., SVM linear and radial at BL -----------------
print("################## FT. SELECTION COMPARISON ANALYSIS (BL) across other models ##################")
type_model_auc_models <- data.frame(
  ftsel = c(rep("DEGs10",4),rep("DEGs100",4), rep("DEGs1000",4), rep("Lasso",4)),
  model = c(rep(c("RandomForest", "LogisticRegression", "SVM-Linear", "SVM-Radial"), 4)),
  AUC = c(c(0.577946127946128, 0.56104797979798, 0.578808922558923,0.567792508417508),
          c(0.585406145,	0.589962121,	0.595875421,	0.594097222),
          c(0.555376684,	0.619307660,	0.631144781,	0.572011785),
          c(0.6310188280,	0.6153631907,	0.6008912656,	0.6401459447))
)
frd <- friedman(groups=type_model_auc_models$ftsel, 
                blocks=type_model_auc_models$model,
                score=type_model_auc_models$AUC)
if (frd$p.value < 0.05) {
  fr_nx1 <- friedman_holm_nx1(groups=type_model_auc_models$ftsel, 
                              blocks=type_model_auc_models$model,
                              score=type_model_auc_models$AUC,
                              col_control=4)
  print(paste("Friedman + Holm post-hoc test: ", names(which.min(fr_nx1$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(fr_nx1$test[2,][which(fr_nx1$test[2,] < 0.05)]), collapse=", ")))
} else{
  print("Friedman hypothesis test was not rejected")
}



# Pairwise comparison: models x models (BL) ------------------------------------
print("################## MODELS COMPARISON ANALYSIS AT BL ##################")
type_model_auc_i <- auc %>%
  dplyr::filter(ftsel=="lasso")
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
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_models_BL.pdf"), width = 14, height = 8)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
}



##### PROBLEM: ONLY 9 GROUPS CAN BE COMPARED IN PAIRWISE COMPARISON
# Pairwise comparison: data types x data types (BL) ----------------------------
print("################## FEATURE TYPES COMPARISON ANALYSIS (BL) ##################")
type_model_auc_i <- auc %>%
  dplyr::filter(ftsel=="lasso")
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
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_datatypes_BL.pdf"), width = 6, height = 4)
  print(fr_nxn$plot)
  dev.off()
  
} else{
  print("Friedman hypothesis test was not rejected")
}


pv.adj <- friedman_holm_nx1(groups=type_model_auc_i$data_type,
                  blocks=type_model_auc_i$model,
                  score=type_model_auc_i$AUC,
                  col_control=7)
# plot 1xn GOBP_mean against all rest---------
df <- reshape2::melt(pv.adj$test[2,])
l <- colnames(pv.adj$test)[order(pv.adj$ranks, decreasing = TRUE)]
df$variable <- factor(df$variable, levels = l)
df <- df[order(df$variable), ]
df$Y <- df$variable[which(is.na(df$value))]
df <- df %>%
  filter(variable !=df$variable[which(is.na(df$value))])
ggplot2::ggplot(df, ggplot2::aes(x = variable, y = Y, 
                                 fill = value)) + ggplot2::geom_tile(col = "white") + 
  ggplot2::scale_fill_continuous("value") + ggplot2::labs(x = "Algorithm", 
                                                            y = "Algorithm") + 
  ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 5, col = "white") + 
  labs(title="Friedman test with adjusted PV (Holm)") + 
  scale_fill_gradientn("Adj. PV" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
  labs(x = "Feature type", y = "Feature type") + theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),    # Set x-axis label size
    axis.title.y = element_text(size = 14),     # Set y-axis label size
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 


# Analyses at BL-TS ------------------------------------------------------------

if (file.exists(IN_FILE_BLTS)) {
  auc_BLTS <- vroom(IN_FILE_BLTS, show_col_types = FALSE)
} else {
  # Data transformations -------------------------------------------------------
  data <- list()
  data[["bl"]] <- read_excel(ORIG_FILE_LASSO, sheet = 8)
  data[["lm_time"]] <- read_excel(ORIG_FILE_LASSO, sheet = 9)
  data[["sd"]] <- read_excel(ORIG_FILE_LASSO, sheet = 10)
  data[["lm_lag"]] <- read_excel(ORIG_FILE_LASSO, sheet = 11)
  
  for (ft in names(data)) {
    data[[ft]] <- data[[ft]] %>%
      fill(1) %>%
      rename_at(1,~"data_type") %>%
      drop_na(2) %>%
      unite(data_type, data_type, 2, sep = "_", remove = TRUE)  %>%
      mutate(data_type = gsub("GENE_-", "GENE", data_type)) %>%
      dplyr::select(-2) %>%
      mutate(across(-1, as.double))
  }
  
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
  
  # export AUC scores table ----------------------------------------------------
  #  merge all AUC scores into full AUC table
  auc_BLTS <- bind_rows(type_model_auc, .id = "fttype")
  # print full AUC table
  readr::write_tsv(auc_BLTS, file = file.path(OUT_DIR, "summary_AUC_scores_BLTS.tsv"))
  
}


# BL vs TS for genes -----------------------------------------------------------
print("################## BL VS TS IN GENES COMPARISON ANALYSIS ##################")

type_model_auc_genes <- auc_BLTS %>%
  dplyr::filter(data_type=="GENE") %>%
  mutate(fttype = toupper(fttype))

print(paste("feature types compared:", paste(unique(type_model_auc_genes$fttype), collapse=", ")))

frd <- friedman(groups=type_model_auc_genes$fttype, 
                blocks=type_model_auc_genes$model,
                score=type_model_auc_genes$AUC)
if (frd$p.value < 0.05) {
  # LM (time) against the other features.
  fr_1xn <- friedman_holm_nx1(groups=type_model_auc_genes$fttype,  
                                        blocks=type_model_auc_genes$model,
                                        score=type_model_auc_genes$AUC, col_control=2)
  
#   fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_genes$fttype,  
#                                        blocks=type_model_auc_genes$model,
#                                        score=type_model_auc_genes$AUC)
#   print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
#               "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
#   # print plot
#   fr_nxn$plot +  
#     scale_x_discrete(labels = function(labels) toupper(labels)) +
#     scale_y_continuous(labels = function(labels) toupper(labels))
#   pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_genes.pdf"), width = 7, height = 4)
#   print(fr_nxn$plot)
#   dev.off()
} else{
  print("Friedman hypothesis test was not rejected")

}



# BL vs TS across data types (best model/data type) ----------------------------
print("################## BL VS TS across data types (best model/data type) COMPARISON ANALYSIS ##################")

type_model_auc_blts <- auc_BLTS %>%
  group_by(data_type, fttype) %>%
  slice_max(AUC) %>%
  mutate(fttype = toupper(fttype))

print(paste("feature types compared:", paste(unique(type_model_auc_blts$fttype), collapse=", ")))

if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts$fttype,  
                                        blocks=type_model_auc_blts$data_type,
                                        score=type_model_auc_blts$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  fr_nxn$plot
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_datatypes.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
  
}

###
# Data types (best model/data type) across BL/TS featurs: NOT POSSIBLE BC > 9 PAIRWISE COMPARISONS.
###
# MEAN, SD, PCA across bl and temporal features (best model/data type)  ------------
print("################## data types across bl and temporal features (best model/data type) COMPARISON ANALYSIS ##################")
type_model_auc_blts_meansd <- type_model_auc_blts %>%
  dplyr::filter(data_type != "GENE") %>%
  dplyr::filter(endsWith(data_type, "mean") | endsWith(data_type, "sd") | endsWith(data_type, "pca"))

print(paste("feature types compared:", paste(unique(type_model_auc_blts_meansd$data_type), collapse=", ")))

frd <- friedman(groups=type_model_auc_blts_meansd$data_type, 
                blocks=type_model_auc_blts_meansd$fttype,
                score=type_model_auc_blts_meansd$AUC)
if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts_meansd$data_type,  
                                        blocks=type_model_auc_blts_meansd$fttype,
                                        score=type_model_auc_blts_meansd$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  fr_nxn$plot +  
    scale_x_discrete(labels = function(labels) toupper(labels)) +
    scale_y_continuous(labels = function(labels) toupper(labels))
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_datatypes.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
  
}



# WITHIN DB: compare aggregations across bl and temporal features (best model/data type)  ------------
for (db in c("GOBP", "GOCC", "CORUM")) {
  print(paste("---------", db, "---------"))
  type_model_auc_blts_db <- type_model_auc_blts %>%
    dplyr::filter(startsWith(data_type, db))
  
  frd <- friedman(groups=type_model_auc_blts_db$data_type, 
                  blocks=type_model_auc_blts_db$fttype,
                  score=type_model_auc_blts_db$AUC)
  
  
  
  if (frd$p.value < 0.05) {
    fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts_db$data_type,  
                                          blocks=type_model_auc_blts_db$fttype,
                                          score=type_model_auc_blts_db$AUC)
    print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
                "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
    # print plot
    fr_nxn$plot +  
      scale_x_discrete(labels = function(labels) toupper(labels)) +
      scale_y_continuous(labels = function(labels) toupper(labels))
    pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_nxn_aggregators_",db,".pdf")), width = 9, height = 4)
    print(fr_nxn$plot)
    dev.off()
  } else{
    print("Friedman hypothesis test was not rejected")
    
  }
  
}


# WITHIN AGGs: compare DB across bl and temporal features (best model/data type)  ------------
for (st in c("mean", "median", "sd", "pathifier", "pca")) {
  print(paste("---------", st, "---------"))
  type_model_auc_blts_st <- type_model_auc_blts %>%
    dplyr::filter(endsWith(data_type, st))
  
  frd <- friedman(groups=type_model_auc_blts_st$data_type, 
                  blocks=type_model_auc_blts_st$fttype,
                  score=type_model_auc_blts_st$AUC)
  
  if (frd$p.value < 0.05) {
    fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts_st$data_type,  
                                          blocks=type_model_auc_blts_st$fttype,
                                          score=type_model_auc_blts_st$AUC)
    print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
                "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
    # print plot
    fr_nxn$plot +  
      scale_x_discrete(labels = function(labels) toupper(labels)) +
      scale_y_continuous(labels = function(labels) toupper(labels))
    pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_nxn_aggregators_",st,".pdf")), width = 9, height = 4)
    print(fr_nxn$plot)
    dev.off()
  } else{
    print("Friedman hypothesis test was not rejected")
    
  }
  
}



# within feature type (e.g. BL), comparison of DB-st RANKS ---------------------


frd <- friedman(groups=type_model_auc_blts$data_type, 
                blocks=type_model_auc_blts$fttype,
                score=type_model_auc_blts$AUC)
RM_TS <- list()
for (tf in c("BL", "LM_TIME", "LM_LAG", "SD")) {
  type_model_auc_blts_r <- type_model_auc_blts %>%
    dplyr::filter(fttype==tf) %>%
    dplyr::select(-model) %>%
    pivot_wider(names_from=data_type, values_from = AUC)
  rm <- colMeans(scmamp::rankMatrix(type_model_auc_blts_r[,-1]))
  rm <- data.frame(rm)
  colnames(rm) <- tf
  RM_TS[[tf]] <- rm
}
rank_df <- do.call(cbind, RM_TS)
rank_df$AVG_RANK <- rowMeans(rank_df)
rank_df <- rank_df[order(rank_df$AVG_RANK), ]
  

# plot ranks
#df <- reshape2::melt(fr_1xn$test[2,])
df <- reshape2::melt(rank_df %>%
  rownames_to_column("data_type") %>%
  pivot_longer(cols= c("BL", "LM_TIME", "LM_LAG", "SD", "AVG_RANK"))) %>%
  select(-variable)
colnames(df) <- c("data_type", "fttype", "value")

l <- rownames(rank_df)[order(rank_df$AVG_RANK, decreasing = TRUE)]
df$data_type <- factor(df$data_type, levels = l)
df <- df[order(df$data_type), ]

rank_plot <- ggplot2::ggplot(df, ggplot2::aes(x = factor(fttype, levels = c("AVG_RANK", "BL", "LM_TIME", "SD", "LM_LAG")), y = data_type, 
                                 fill = value)) + ggplot2::geom_tile(col = "white") + 
  ggplot2::scale_fill_continuous("value") + ggplot2::labs(x = "Algorithm", 
                                                          y = "Algorithm") + 
  ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 5, col = "white") + 
  labs(title="Ranked performance of database-pooling operator combinations across baseline and temporal features") + 
  scale_fill_gradientn("Rank" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
  labs(x = "Feature type", y = "Data type") + theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),    # Set x-axis label size
    axis.title.y = element_text(size = 14),     # Set y-axis label size
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 

# print plot
pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_ranks_aggs_BLTS.pdf")), width = 12, height = 7)
print(rank_plot)
dev.off()



# plot 1xn GOBP_mean against all rest---------
fr_1xn <- friedman_holm_nx1(groups=type_model_auc_blts$data_type,
                            blocks=type_model_auc_blts$fttype,
                            score=type_model_auc_blts$AUC,
                            col_control=7)

df <- reshape2::melt(fr_1xn$test[2,])
l <- colnames(fr_1xn$test)[order(fr_1xn$ranks, decreasing = TRUE)]
df$variable <- factor(df$variable, levels = l)
df <- df[order(df$variable), ]
df$Y <- df$variable[which(is.na(df$value))]
df <- df %>%
  filter(variable !=df$variable[which(is.na(df$value))])
plot <- ggplot2::ggplot(df, ggplot2::aes(x = variable, y = Y, 
                                 fill = value)) + ggplot2::geom_tile(col = "white") + 
  ggplot2::scale_fill_continuous("value") + ggplot2::labs(x = "Algorithm", 
                                                          y = "Algorithm") + 
  ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 5, col = "white") + 
  labs(title="") + 
  scale_fill_gradientn("Adj. PV" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
  labs(x = "Data type", y = "") + theme(
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


# print plot
pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_1xn_GOBP_mean.pdf")), width = 14, height = 2)
print(plot)
dev.off()





 # within feature type (e.g. BL), comparison of DB-st using mean and sd?
# for the same aggregations: compare DB.
type_model_auc_blts_ts <- type_model_auc_blts %>%
  dplyr::filter(data_type != "GENE") %>%
  dplyr::filter(fttype == "SD") %>%
  separate(data_type, into = c("db", "st"), sep = "_")
  

print(paste("feature types compared:", paste(unique(type_model_auc_blts_ts$st), collapse=", ")))

frd <- friedman(groups=type_model_auc_blts_ts$db, 
                blocks=type_model_auc_blts_ts$st,
                score=type_model_auc_blts_ts$AUC)
if (frd$p.value < 0.05) {
  fr_1xn <- friedman_holm_nx1(groups=type_model_auc_blts$db,  
                              blocks=type_model_auc_blts$st,
                              score=type_model_auc_blts$AUC, col_control=7)
  
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts_ts$db,  
                                        blocks=type_model_auc_blts_ts$st,
                                        score=type_model_auc_blts_ts$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  fr_nxn$plot +  
    scale_x_discrete(labels = function(labels) toupper(labels)) +
    scale_y_continuous(labels = function(labels) toupper(labels))
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_datatypes.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
  
}


