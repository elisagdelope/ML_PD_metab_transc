# Title: lx_nCV_rf_MLclass.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs nested crossvalidation for a random forest classifier on RNAseq static data (BL).
# Usage (gene level): Rscript lx_nCV_rf_MLclass.R 
# Usage (aggregated level): Rscript lx_nCV_rf_MLclass.R -l pw -s mean 
# Data: data from metabolomics at specific timepoint (V0), clinical data (diagnosis).

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
library(fastshap)
library(randomForest)


# 1 time analysis: to do on BL and TS data (correlation and near-zero filter already applied to the input data).

# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-V0-PD"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_PLOTS <- paste0("../data/", analysis_name , "/01-plots") 
OUT_DIR_PLOTS_SHAP <- paste0(OUT_DIR_PLOTS, "/shap_data") 
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv")
IN_DIR <- "../data/00-cleansing/"
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")
target = "DIAGNOSIS" # pass as cmd parameter
myseed = 111

# Add command line arguments
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (pw / metab)", default = "metab", type = "string",  nargs='?')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string",  nargs='?')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene / aggregations
st = tolower(argv$stat) # stat at aggregation level


if (e_level == "METAB") { # normalized counts at gene level
  metab_4ML.FILE <- file.path(OUT_DIR, paste0("data_cv_metab_4ML_", target, ".tsv"))
  hout_4ML.FILE  <- file.path(OUT_DIR, paste0("data_test_metab_4ML_", target, ".tsv"))
  OUT.AUC.FILE <- file.path(OUT_DIR, paste0("lt_", target, "_rf_AUC.tsv"))
  OUT.ERROR.FILE <- file.path(OUT_DIR, paste0("lt_", target, "_rf_ERROR.tsv"))
  OUT.MODEL.FILE <- file.path(OUT_DIR, paste0("lt_", target, "_rf_final_model.rds"))
  OUT.PROBS.FILE <- file.path(OUT_DIR, paste0("lt_", target, "_rf_final_probs.tsv"))
  OUT.FPERF.FILE <- file.path(OUT_DIR, paste0("lt_", target, "_rf_final_performance.tsv"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste0("lt_", target, "_rf_varImp.pdf"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", target, "_rf_shap_values.csv"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", target, "_rf_shap_data.csv"))
  features_varname <- "METABOLITES_ID"
  
} else { # aggregations
  metab_4ML.FILE <- file.path(OUT_DIR_PATHWAY, paste0("PW_", st, "_data_cv_metab_4ML_", target, ".tsv"))
  hout_4ML.FILE  <- file.path(OUT_DIR_PATHWAY, paste0("PW_", st, "_data_test_metab_4ML_", target, ".tsv"))
  OUT.AUC.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", target, "_rf_AUC.tsv"))
  OUT.ERROR.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", target, "_rf_ERROR.tsv"))
  OUT.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", target, "_rf_final_model.rds"))
  OUT.PROBS.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", target, "_rf_final_probs.tsv"))
  OUT.FPERF.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", target, "_rf_final_performance.tsv"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste0("lt_PW_", st, "_", target, "_rf_varImp.pdf"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", target, "_rf_shap_values.csv"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", target, "_rf_shap_data.csv"))
  features_varname <- "PATHWAY_NAME"
  
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(metab_4ML.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
  }
}


# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY)) | (!dir.exists(OUT_DIR_PLOTS)) | (!dir.exists(OUT_DIR_PLOTS_SHAP))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
  dir.create(OUT_DIR_PLOTS_SHAP, recursive = T)
}



# Data load --------------------------------------------------------------------
# log-transformed peak area & clinical data.
annotation <- vroom(ANNOTATION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 
metab_4ML <- vroom(metab_4ML.FILE, col_types = cols()) 
hout_4ML <- vroom(hout_4ML.FILE, col_types = cols()) 



# Data reformatting for ML -----------------------------------------------------

# kNN imputation for BMI variable
pheno <- VIM::kNN(pheno, variable = "BMI", k= 5, imp_var = F)
pheno <- pheno %>%
  filter(SAMPLE_ID %in% metab_4ML$SAMPLE_ID) %>% # filter samples in training cv data
  mutate_at(target, factor)

# flip the dataset as rows should correspond to metabolites and columns to samples
metab.t = metab_4ML %>% 
  dplyr::select(-any_of(target)) %>%
  pivot_longer(!SAMPLE_ID, names_to = features_varname, values_to = "COUNT") %>% 
  pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
  column_to_rownames(features_varname)

metab_4ML <- metab_4ML %>%
  column_to_rownames(var = "SAMPLE_ID")
hout_4ML <- hout_4ML %>%
  column_to_rownames(var = "SAMPLE_ID")



# create nested cv -------------------------------------------------------------
set.seed(myseed)
split_index_outer <- createDataPartition(metab_4ML[[target]], p = 0.75, list = FALSE, times = 10)
head(split_index_outer, 10)

# placeholders for performance measures
ERROR.OUTER <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_outer)))
AUC.OUTER <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_outer)))
mymodel <- "rf"
colnames(ERROR.OUTER) <- c('fold', mymodel)
colnames(AUC.OUTER) <- c('fold', mymodel)

for (i in 1:ncol(split_index_outer)){
  
  # use ith column of split_index to create feature and target training/test sets (subjects split)
  metab.t_train <- metab.t[, split_index_outer[,i]] 
  pheno_train <- pheno[split_index_outer[,i], ]
  
  XY_train <- metab_4ML[split_index_outer[,i],] %>%  # keeps the target var 
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(metab_4ML[[target]])))))) # target values to 'valid' names: {X0, X1}
  X_test  <- metab_4ML[-split_index_outer[,i], !(names(metab_4ML) %in% c(target))]  # without target var
  Y_test <- metab_4ML[-split_index_outer[,i], ] %>% 
    dplyr::select(all_of(target)) %>%  
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(XY_train[[target]])))))) # target values to 'valid' names: {X0, X1}
  
  print(paste("Dimensions of outer train set", as.character(dim(XY_train)[1]), "x", as.character(dim(XY_train)[2]-1), sep = " "))
  print(paste("Dimensions of outer test set", as.character(dim(X_test)[1]), "x", as.character(dim(X_test)[2]), sep = " "))
  
  ## no downsampling needed
  
  fitControl <- trainControl(method = "repeatedcv",  # inner cv loop 2
                             number = 5, # 5-fold cv
                             repeats = 2,
                             p = 0.7,
                             verboseIter = FALSE,
                             search = "grid",
                             savePredictions = "final",
                             classProbs = TRUE,
                             # allowParallel = TRUE,
                             summaryFunction = twoClassSummary) 
  
  model_fit <- caret::train(form = as.formula(paste(target, "~ .")), 
                                              data = XY_train, 
                                              method = mymodel, 
                                              trControl = fitControl,
                                              metric = "ROC",
                                              verbose = FALSE,
                                              preProcess = c('scale', 'center')) # scale & center data
  
  pred_test = predict(model_fit, newdata = X_test, type = "prob") 

  # track outer cv performance: mean error, AUC
  error_pred <- as.factor(predict(model_fit, newdata = X_test))
  print(caret::confusionMatrix(error_pred, Y_test %>% pull(), positive='X1'))
  error_pred <- mean(ifelse(Y_test %>% pull() != error_pred, 1, 0))
  ERROR.OUTER[i, mymodel] <- error_pred
  ERROR.OUTER[i, 'fold'] <- i
  
  auc_pred <- pROC::roc(response = as.vector(as.matrix(Y_test)), predictor = pred_test[[2]], directon = "<")$auc
  AUC.OUTER[i, mymodel] <- auc_pred
  AUC.OUTER[i, 'fold'] <- i
}

# Add {avg, sd} performance of each model in the outer sets (across outer folds)
ERROR.OUTER <- rbind(ERROR.OUTER, c("avg", mean(ERROR.OUTER[1:ncol(split_index_outer), 2])))
ERROR.OUTER[, 2] <- as.numeric(ERROR.OUTER[, 2] )
ERROR.OUTER <- rbind(ERROR.OUTER, c("sd", colSds(as.matrix(ERROR.OUTER[1:ncol(split_index_outer), 2]))))
ERROR.OUTER[, 2] <- as.numeric(ERROR.OUTER[, 2] )
AUC.OUTER <- rbind(AUC.OUTER, c("avg", mean(AUC.OUTER[1:ncol(split_index_outer), 2])))
AUC.OUTER[, 2] <- as.numeric(AUC.OUTER[, 2] )
AUC.OUTER <- rbind(AUC.OUTER, c("sd", colSds(as.matrix(AUC.OUTER[1:ncol(split_index_outer), 2]))))
AUC.OUTER[, 2] <- as.numeric(AUC.OUTER[, 2] )

print(knitr::kable(ERROR.OUTER, row.names = F))
print(knitr::kable(AUC.OUTER, row.names = F))

# export outputs  --------------------------------------------------------------
readr::write_tsv(ERROR.OUTER, file = OUT.ERROR.FILE) 
readr::write_tsv(AUC.OUTER, file = OUT.AUC.FILE) 
rm(XY_train, X_test, Y_test)





####################
# BUILD FINAL MODEL
####################

metab_4ML <- metab_4ML %>%
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(metab_4ML[[target]]))))))

## downsampling not needed

# train final model
fitControl <- trainControl(method = "repeatedcv",  
                           number = 10, # 10-fold cv
                           repeats = 2,
                           p = 0.7,
                           verboseIter = FALSE,
                           search = "grid",
                           savePredictions = "final",
                           classProbs = TRUE,
                           # allowParallel = TRUE,
                           summaryFunction = twoClassSummary) 

f_model <- caret::train(form = as.formula(paste(target, "~ .")), 
                        data = metab_4ML, 
                        method = mymodel, 
                        trControl = fitControl,
                        metric = "ROC",
                        verbose = FALSE,
                        preProcess = c('scale', 'center')) # scale & center data

# final assessment on held-out data
X_hout <- hout_4ML[, !(names(hout_4ML) %in% c(target))]  
Y_hout <- hout_4ML %>% 
  dplyr::select(all_of(target)) %>% 
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(hout_4ML[[target]]))))))


Y_hat = predict(f_model, newdata = X_hout, type = "prob")
colnames(Y_hat) <- paste(mymodel, colnames(Y_hat), sep = ".")
colnames(Y_hat) <- sub("X", "", colnames(Y_hat))
rownames(Y_hat) <- rownames(X_hout)

error_pred <- as.factor(predict(f_model, newdata = X_hout))
print(caret::confusionMatrix(error_pred, Y_hout %>% pull(), positive='X1'))
error_pred <- mean(ifelse(Y_hout %>% pull() != error_pred, 1, 0))
auc_pred <- pROC::roc(response = as.vector(as.matrix(Y_hout)), predictor = Y_hat[[2]])$auc

info_df <- data.frame("AUC" = auc_pred,
                      "Error" = error_pred)
print(info_df)



# export outputs  --------------------------------------------------------------
saveRDS(f_model, OUT.MODEL.FILE) # save the model object
readr::write_tsv(Y_hat %>% rownames_to_column("PATIENT_ID"), file = OUT.PROBS.FILE) 
readr::write_tsv(info_df, file = OUT.FPERF.FILE) 

# plots
pdf(file = OUT.PLOTS.FILE, width = 14, height = 8)
# if model has its own variable scoring:
variable_imp <- varImp(f_model, scale = FALSE) # area under the ROC curve is computed for each class
if (("Overall" %in% names(variable_imp$importance)) & (nrow(variable_imp$importance) > 20)) {
  vi <- as.data.frame(variable_imp$importance) %>%
    rownames_to_column("PATHWAY_NAME") %>%
    arrange(desc(Overall)) %>%
    top_n(20) %>%
    mutate_at(vars(PATHWAY_NAME), 
              ~gsub("`", "", .)) %>%
    column_to_rownames("PATHWAY_NAME")
  variable_imp$importance <-vi
  plot(variable_imp, main=paste(mymodel, "- Variable Importance"))
}
variable_imp_r <- caret::varImp(f_model, useModel = FALSE, scale = FALSE) #from caret 
if (nrow(variable_imp_r$importance) > 20) {
  vi <- as.data.frame(variable_imp_r$importance) %>%
    rownames_to_column("PATHWAY_NAME") %>%
    arrange(desc(X0)) %>%
    top_n(20) %>%
    mutate_at(vars(PATHWAY_NAME), 
              ~gsub("`", "", .)) %>%
    column_to_rownames("PATHWAY_NAME")
  variable_imp_r$importance <- vi
}
plot(variable_imp_r, main=paste(mymodel, "- Variable Importance (ROC analysis)"))
dev.off()

readr::write_csv(metab_4ML[,-ncol(metab_4ML)], file = OUT.SHAPDATA.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

