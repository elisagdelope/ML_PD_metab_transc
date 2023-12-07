# Title: ppmi_nCV_rf_MLclass.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs nested crossvalidation for a random forest classifier on RNAseq static data (BL).
# Usage (gene level): Rscript ppmi_nCV_rf_MLclass.R filtercor (for applying correlation filter)
# Usage (aggregated level): Rscript ppmi_nCV_rf_MLclass.R filtercor -l gobp -s mean 
# Data: data from expression counts at specific timepoint (BL), clinical data (diagnosis).

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


# 1 time analysis: to do on BL and TS data with and without correlation filter.
# test the variable importance plots -> top 50 ? 30? 20? features?

# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-BL-PD"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_PLOTS <- paste0("../data/", analysis_name , "/01-plots") 
OUT_DIR_PLOTS_SHAP <- paste0(OUT_DIR_PLOTS, "/shap_data") 
PHENO.FILE <- file.path(OUT_DIR, "pheno_BL.tsv")
target = "DIAGNOSIS" # pass as cmd parameter
myseed = 111

# Add command line arguments
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "filtercor", help = "boolean string (yes/no) whether to use a correlation filter (pearson, 0.85)", flag=TRUE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string",  nargs='?')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string",  nargs='?')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
flt_cor = toupper(argv$filtercor) # trigger for corr filter
e_level = toupper(argv$level) # gene / aggregations
st = tolower(argv$stat) # stat at aggregation level

if (flt_cor == "FALSE") {
  if (e_level == "GENE") { # normalized counts at gene level
    EXPRS.FILE <- file.path(OUT_DIR, "flt_norm_star_BL.tsv")
    OUT.AUC.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_AUC.tsv", sep = "_"))
    OUT.ERROR.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_ERROR.tsv", sep = "_"))
    OUT.MODEL.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_final_model.rds", sep = "_"))
    OUT.PROBS.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_final_probs.tsv", sep = "_"))
    OUT.FPERF.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_final_performance.tsv", sep = "_"))
    OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, "norm", target, "rf_varImp.pdf", sep = "_"))
    OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", target, "rf_shap_values.csv", sep = "_"))
    OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", target, "rf_shap_data.csv", sep = "_"))
    features_varname <- "GENEID"
    
  } else { # aggregations
    EXPRS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
    OUT.AUC.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_AUC.tsv", sep = "_"))
    OUT.ERROR.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_ERROR.tsv", sep = "_"))
    OUT.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_final_model.rds", sep = "_"))
    OUT.PROBS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_final_probs.tsv", sep = "_"))
    OUT.FPERF.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_final_performance.tsv", sep = "_"))
    OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, st, target, "rf_varImp.pdf", sep = "_"))
    OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, target, "rf_shap_values.csv", sep = "_"))
    OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, target, "rf_shap_data.csv", sep = "_"))
    features_varname <- paste0(e_level, "_NAME")
    
    if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier")) | (!file.exists(EXPRS.FILE))) { 
      stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
    }
  }
} else {
  if (e_level == "GENE") { # normalized counts at gene level
    EXPRS.FILE <- file.path(OUT_DIR, "flt_norm_star_BL.tsv")
    OUT.AUC.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_fltcor_AUC.tsv", sep = "_"))
    OUT.ERROR.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_fltcor_ERROR.tsv", sep = "_"))
    OUT.MODEL.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_fltcor_final_model.rds", sep = "_"))
    OUT.PROBS.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_fltcor_final_probs.tsv", sep = "_"))
    OUT.FPERF.FILE <- file.path(OUT_DIR, paste(e_level, "norm", target, "rf_fltcor_final_performance.tsv", sep = "_"))
    OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, "norm", target, "rf_fltcor_varImp.pdf", sep = "_"))
    OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", target, "rf_fltcor_shap_values.csv", sep = "_"))
    OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", target, "rf_fltcor_shap_data.csv", sep = "_"))
    features_varname <- "GENEID"
    
  } else { # aggregations
    EXPRS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
    OUT.AUC.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_fltcor_AUC.tsv", sep = "_"))
    OUT.ERROR.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_fltcor_ERROR.tsv", sep = "_"))
    OUT.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_fltcor_final_model.rds", sep = "_"))
    OUT.PROBS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_fltcor_final_probs.tsv", sep = "_"))
    OUT.FPERF.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, target, "rf_fltcor_final_performance.tsv", sep = "_"))
    OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, st, target, "rf_fltcor_varImp.pdf", sep = "_"))
    OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, target, "rf_fltcor_shap_values.csv", sep = "_"))
    OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, target, "rf_fltcor_shap_data.csv", sep = "_"))
    features_varname <- paste0(e_level, "_NAME")
    
    if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier")) | (!file.exists(EXPRS.FILE))) { 
      stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
    }
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
# star raw counts & (pre-filtered) clinical data.
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddlf")) 
expr <- vroom(EXPRS.FILE, col_types = cols()) %>%
  rename_with(toupper)
expr = expr[!duplicated(expr[[features_varname]]),]
##### TO BE REMOVED #####
#pheno <- pheno[1:200,] # 200 subjects, 1000 genes
#expr <- expr[,1:201]
#expr <- expr[1:500,]
#########################

# Data reformatting for ML -----------------------------------------------------

expr_4ML <- expr %>%
  pivot_longer(cols = stringr::str_subset(colnames(expr), "[0-9]{4}\\.[A-Z]"), names_to = "SAMPLE", values_to = "EXPRS") %>%
  separate(SAMPLE, c("PATIENT_ID", "VISIT"), sep = '\\.') %>%
  pivot_wider(names_from = all_of(features_varname), values_from = EXPRS) %>%
  dplyr::select(-VISIT)

pheno_4ML <- pheno %>%
  dplyr::select(all_of(c(target, 'PATIENT_ID'))) %>%
  mutate(!!target := case_when(get(target) == "HC" ~ 0,
                               get(target) == "PD" ~ 1))

expr_4ML <- expr_4ML %>%
  inner_join(pheno_4ML, 
             by = "PATIENT_ID") %>%
  mutate_at(target, factor)
rm(expr, pheno, pheno_4ML)



# apply unsupervised filters ---------------------------------------------------
n <- ncol(expr_4ML)
# remove zero variance features 
nzv <- nearZeroVar(expr_4ML[, !(names(expr_4ML) %in% c(target))], freqCut = 10 ) 
if (length(nzv) != 0 ) {
  expr_4ML <- expr_4ML[,-nzv]
}

if (flt_cor == "TRUE") {
  # remove highly correlated features 
  cor_df = cor(expr_4ML[,-c(1,ncol(expr_4ML))]) # remove patient ID and diagnosis variables
  hc = findCorrelation(cor_df, cutoff=0.85, names = TRUE) 
  hc = sort(hc)
  expr_4ML = expr_4ML[,-which(names(expr_4ML) %in% c(hc))]
  print("NZV and correlation filters successfully applied")
  print(paste((length(hc) + length(nzv)), "features were removed out of", n))
} else {
  print("NZV filter successfully applied")
  print(paste(length(nzv), "features were removed out of", n))
}

expr_4ML <- expr_4ML %>% 
  column_to_rownames("PATIENT_ID")



# create training/held-out set -------------------------------------------------
set.seed(myseed-1)
inTraining <- createDataPartition(expr_4ML[[as.character(target)]], p = .85, list = FALSE)
hout_4ML  <- expr_4ML[-inTraining, ]
expr_4ML <- expr_4ML[inTraining, ] 


# create nested cv -------------------------------------------------------------
set.seed(myseed)
split_index_outer <- createDataPartition(expr_4ML[[target]], p = 0.7, list = FALSE, times = 10)
head(split_index_outer, 10)

# placeholders for performance measures
ERROR.OUTER <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_outer)))
AUC.OUTER <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_outer)))
mymodel <- "rf"
colnames(ERROR.OUTER) <- c('fold', mymodel)
colnames(AUC.OUTER) <- c('fold', mymodel)

for (i in 1:ncol(split_index_outer)){
  
  # use ith column of split_index to create training/test sets (subjects split)
  XY_train <- expr_4ML[split_index_outer[,i],] %>%  # keeps the target var 
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(expr_4ML[[target]])))))) # target values to 'valid' names: {X0, X1}
  
  X_test  <- expr_4ML[-split_index_outer[,i], !(names(expr_4ML) %in% c(target))]  # without target var
  Y_test <- expr_4ML[-split_index_outer[,i], ] %>% 
    dplyr::select(all_of(target)) %>%  
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(expr_4ML[[target]])))))) # target values to 'valid' names: {X0, X1}
  
  print(paste("Dimensions of outer train set", as.character(dim(XY_train)[1]), "x", as.character(dim(XY_train)[2]-1), sep = " "))
  print(paste("Dimensions of outer test set", as.character(dim(X_test)[1]), "x", as.character(dim(X_test)[2]), sep = " "))
  
  # downsampling 
  set.seed(myseed + 1)
  XY_train <- downSample(x=XY_train[,-ncol(XY_train)] %>% rownames_to_column(var = "PATIENT_ID"),
                         y=XY_train[[target]],
                         yname = target) %>%
    arrange(PATIENT_ID) %>% 
    column_to_rownames(var = "PATIENT_ID")
  print("downsampling in outer loop successfully applied:")
  print(table(XY_train[[target]]))
  
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
                                              verbosity = 0) # parameter for xgb-trees models
  
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

# downsampling 
set.seed(myseed + 3)
expr_4ML <- downSample(x=expr_4ML[,-ncol(expr_4ML)] %>% rownames_to_column(var = "PATIENT_ID"),
                       y=expr_4ML[[target]],
                       yname = target) %>%
  arrange(PATIENT_ID) %>% 
  column_to_rownames(var = "PATIENT_ID") %>%
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(expr_4ML[[target]]))))))

print("downsampling in final training successfully applied:")
print(table(expr_4ML[[target]]))

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
                        data = expr_4ML, 
                        method = mymodel, 
                        trControl = fitControl,
                        metric = "ROC",
                        verbose = FALSE,
                        verbosity = 0) # parameter for xgb-trees models

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
variable_imp = varImp(f_model, scale = FALSE) # area under the ROC curve is computed for each class
# sort variable_imp$importance only the top 20
variable_imp$importance <- variable_imp$importance %>% 
                            arrange(-Overall) %>%
                            top_n(20)
  
plot(variable_imp, main=paste(mymodel, "- Variable Importance"))
variable_imp_r <- caret::varImp(f_model, useModel = FALSE, scale = FALSE) #from caret 
variable_imp_r$importance <- variable_imp_r$importance  %>% 
                              arrange(-X0) %>%
                              top_n(20)
plot(variable_imp_r, main=paste(mymodel, "- Variable Importance (ROC analysis)"))
dev.off()

# Prediction wrapper
pfun <- function(object, newdata) {
  predict(f_model, newdata = newdata, type = "prob")[,2]
}
shap <- fastshap::explain(f_model, X = expr_4ML[,-ncol(expr_4ML)], pred_wrapper = pfun, nsim = 10)

readr::write_csv(shap, file = OUT.SHAPVALUES.FILE)
readr::write_csv(expr_4ML[,-ncol(expr_4ML)], file = OUT.SHAPDATA.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

