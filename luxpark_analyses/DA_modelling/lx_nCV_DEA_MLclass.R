# Title: lx_nCV_DEA_MLclass.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs nested crossvalidation for feature selection (through DA) + binary ML classifiers on metabolomics static data (V0).
# Usage (gene level): Rscript lx_nCV_DEA_MLclass.R 10 (for 10 features)
# Usage (aggregated level): Rscript lx_nCV_DEA_MLclass.R 10 -l pw -s mean 
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
library(limma)
library(argparser, quietly = TRUE)
library(matrixStats)
library(randomForest)
library(glmnet)
library(LiblineaR)
library(caTools)
library(e1071) 
library(kernlab)
library(xgboost)
library(RSNNS)
library(deepboost)
library(VIM)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-V0-PD" # pass as cmd parameter, typically "02-pred-V0-PD" | "02-pred-V0-UPDRS3-class"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_PLOTS <- paste0("../data/", analysis_name , "/01-plots") 
OUT_DIR_PLOTS_SHAP <- paste0(OUT_DIR_PLOTS, "/shap_data") 
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv") # lx_pheno.tsv or pheno_V0.tsv
IN_DIR <- "../data/00-cleansing/"
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")
M1342.FILE <- file.path(IN_DIR, "M1342.tsv") # M1342 is ID for 3-methoxytyrosine: confounder proxy

target = "DIAGNOSIS" # pass as cmd parameter, typically "DIAGNOSIS" | "UPDRS__3_binary"
myseed = 111
source("deepboost_prob.R")

# Add command line arguments
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "k", help = "number of features selected from DEA", default = 10, short = "k")
p <- add_argument(parser = p, arg = "--level", help = "level of features (pw / metab)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
k = as.integer(argv$k) # number of features to select
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level
print(argv)


if (e_level == "METAB") {
  METAB.FILE <- file.path(OUT_DIR, "log_transformed_V0.tsv")
  OUT.AUC.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_AUC.tsv"))
  OUT.ERROR.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_ERROR.tsv"))
  OUT.MODEL.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_final_model.rds"))
  OUT.PROBS.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_final_probs.tsv"))
  OUT.FPERF.FILE <- file.path(OUT_DIR, paste0("lt_", k, "ft_", target, "_final_performance.tsv"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste0("lt_", k, "ft_", target, "_varImp.pdf"))
#  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", k, "ft_", target, "_shap_values.csv"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_", k, "ft_", target, "_shap_data.csv"))
  features_varname <- "METABOLITES_ID"
} else { # aggregations
  METAB.FILE <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, "_V0.tsv"))
  OUT.AUC.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_AUC.tsv"))
  OUT.ERROR.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_ERROR.tsv"))
  OUT.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_final_model.rds"))
  OUT.PROBS.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_final_probs.tsv"))
  OUT.FPERF.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", k, "ft_", target, "_final_performance.tsv"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste0("lt_PW_", st, "_", k, "ft_", target, "_varImp.pdf"))
#  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", k, "ft_", target, "_shap_values.csv"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste0("lt_PW_", st, "_", k, "ft_", target, "_shap_data.csv"))
  features_varname <- "PATHWAY_NAME"
  
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(METAB.FILE))) { 
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
M1342 <- vroom(M1342.FILE, col_types = cols()) 
pheno <- vroom(PHENO.FILE, col_types = cols())
metab <- vroom(METAB.FILE, col_types = cols()) 
metab = metab[!duplicated(metab[["SAMPLE_ID"]]), ]



# Data reformatting for ML -----------------------------------------------------

# kNN imputation for BMI variable
pheno <- VIM::kNN(pheno, variable = "BMI", k= 5, imp_var = F)

if(!all(pheno[[target]] %in% c(0,1))) {
  pheno_4ML <- pheno %>%
    dplyr::select(all_of(c(target, 'SAMPLE_ID'))) %>%
    mutate(!!target := case_when(get(target) == "HC" ~ 0,
                                 get(target) == "PD" ~ 1))
  pheno <- pheno %>% 
    mutate_at(target, factor)
} else {
  pheno_4ML <- pheno %>%
    dplyr::select(all_of(c(target, 'SAMPLE_ID')))
  pheno <- pheno %>% 
    mutate_at(target, factor) %>%
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(pheno[[target]]))))))
  
}

dim_metab <- dim(metab)



# apply unsupervised filters ---------------------------------------------------

# remove near zero variance features 
nzv = nearZeroVar(metab[,-c(1:3)], names = TRUE)
if (length(nzv) > 0) {
  metab <- metab %>%
    dplyr::select(-any_of(nzv))
}

# remove highly correlated features 
cor_df = cor(metab[,-c(1:3)]) # remove patient ID variables
hc = findCorrelation(cor_df, cutoff=0.85, names = TRUE) 
hc = sort(hc)
if (length(hc) > 0) {
  metab = metab[,-which(names(metab) %in% c(hc))]
}

print("NZV, correlation filters and treatment effects successfully applied")
print(paste((length(hc) + length(nzv)), "features were removed out of", (dim_metab[2]-3)))

# flip the dataset as rows should correspond to metabolites and columns to samples
metab.t = metab %>% 
  dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>%
  pivot_longer(!SAMPLE_ID, names_to = features_varname, values_to = "COUNT") %>% 
  pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
  column_to_rownames(features_varname)

# add target variable information
metab_4ML <- metab %>%
  dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>%
  inner_join(pheno_4ML, 
             by = "SAMPLE_ID") %>%
  mutate_at(target, factor) %>% 
  column_to_rownames("SAMPLE_ID")

rm(pheno_4ML, metab)



# create training/held-out set -------------------------------------------------
set.seed(myseed-1)
inTraining <- createDataPartition(metab_4ML[[as.character(target)]], p = .85, list = FALSE)
hout_4ML  <- metab_4ML[-inTraining, ]
metab_4ML <- metab_4ML[inTraining, ] 
metab.t <- metab.t[, inTraining]
pheno <- pheno[inTraining, ]

#### export data 
#if (e_level == "METAB") {
#  readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR, paste0("data_cv_metab_4ML_", target, ".tsv")))
#  readr::write_tsv(hout_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR, paste0("data_test_metab_4ML_", target, ".tsv")))
#} else {
#  readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_PATHWAY, paste0(e_level, "_", st, "_data_cv_metab_4ML_", target, ".tsv")))
#  readr::write_tsv(hout_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_PATHWAY, paste0(e_level, "_", st, "_data_test_metab_4ML_", target, ".tsv")))
#}
####



# create nested cv -------------------------------------------------------------
set.seed(myseed)
split_index_outer <- createDataPartition(metab_4ML[[target]], p = 0.8, list = FALSE, times = 10)
head(split_index_outer, 10)

# define models & parameters
MODELS = list()
model_list_caret <- c("rf", "glmnet", 'regLogistic', "svmLinear2", 
                      "svmRadial", "svmPoly", "xgbDART", 'xgbLinear', "deepboost", "mlpWeightDecayML") # removed 'LogitBoost' as it performed poorly
parameters <- c(seq(0.001, 0.1, by = 0.01), seq(0.1, 2, by = 0.15),  seq(2.5, 5, 0.5))

# placeholders for features and performance measures
FEAT.SELECTION = list()
ERROR <- list()
AUC <- list()
ERROR.OUTER <- data.frame(matrix(ncol = length(model_list_caret) +1, nrow = ncol(split_index_outer)))
AUC.OUTER <- data.frame(matrix(ncol = length(model_list_caret) +1, nrow = ncol(split_index_outer)))
colnames(ERROR.OUTER) <- c('fold', model_list_caret)
colnames(AUC.OUTER) <- c('fold', model_list_caret)

for (i in 1:ncol(split_index_outer)){
  FEAT.SELECTION[[i]] = list()
  MODELS[[i]] <- list()
  
  # use ith column of split_index to create feature and target training/test sets (subjects split)
  metab.t_train <- metab.t[, split_index_outer[,i]] 
  pheno_train <- pheno[split_index_outer[,i], ]
  
  XY_train <- metab_4ML[split_index_outer[,i],] # keeps the target var
  X_test  <- metab_4ML[-split_index_outer[,i], !(names(metab_4ML) %in% c(target))]  # without target var
  Y_test <- metab_4ML[-split_index_outer[,i], ] %>% 
    dplyr::select(all_of(target)) %>%  
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(XY_train[[target]])))))) # target values to 'valid' names: {X0, X1}
  
  print(paste("Dimensions of outer train set", as.character(dim(XY_train)[1]), "x", as.character(dim(XY_train)[2]-1), sep = " "))
  print(paste("Dimensions of outer test set", as.character(dim(X_test)[1]), "x", as.character(dim(X_test)[2]), sep = " "))
  
  ## no downsampling needed
  
  # inner loop partitions
  set.seed(myseed + 2)
  split_index_inner <- createDataPartition(XY_train[[target]], p = 0.8, list = FALSE, times = 5)
  head(split_index_inner, 10)  
  
  
  for (model in model_list_caret) {
    print(paste("Fitting model", model, "..."))
    
    ERROR[[model]] <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_inner)))
    AUC[[model]] <- data.frame(matrix(ncol = 2, nrow = ncol(split_index_inner)))
    colnames(ERROR[[model]]) <- c('fold', model)
    colnames(AUC[[model]]) <- c('fold', model)
    MODELS[[i]][[model]]<- list()
    
    for(j in 1:ncol(split_index_inner)){
      
      if (which(model == model_list_caret) == 1) {
        # DEA to select features -----
        pheno_train_i <- pheno_train[split_index_inner[,j],]
        M1342_i <- M1342[M1342$SAMPLE_ID %in% pheno_train_i$SAMPLE_ID, ]$M1342 # confounding factor 3-methoxytyrosine
        metab.t_train_i <- metab.t_train[, split_index_inner[,j]]
          
        # Apply Bayesian Moderated t-statistic -----------------------------------------
        # empirical bayes t-test
        design <- model.matrix(~0 + pheno_train_i[[target]] + pheno_train_i$AGE + pheno_train_i$GENDER + pheno_train_i$BMI + M1342_i) # diagnosis + confounding factors
        colnames(design) <- c(levels(pheno_train_i[[target]]), "AGE", "GENDER", "BMI", "M1342")
        contrast.matrix <- limma::makeContrasts(paste0(levels(pheno_train_i[[target]])[1], '-', levels(pheno_train_i[[target]])[2]), levels = design) # comparison for diagnosis
        
        xfit <- limma::lmFit(metab.t_train_i, design)
        xfit <- limma::contrasts.fit(xfit, contrast.matrix)
        ebayes <- limma::eBayes(xfit)
        lm_summary <- summary(decideTests(ebayes))
        print(lm_summary)
        lm_res <- topTable(ebayes, coef=1, number = nrow(ebayes), sort.by = "P")
        lm_res <- as.data.frame(lm_res) %>% 
          rownames_to_column(var = features_varname) %>%
          filter(!is.na(adj.P.Val)) %>%
          arrange(P.Value) %>%
          slice_head(n = k)
        
        feature_ids <- lm_res[[1]] 
        FEAT.SELECTION[[i]][[j]] <- lm_res[[1]]
        print(paste(length(feature_ids), "features were selected by DEA"))
      } else {
        feature_ids <- FEAT.SELECTION[[i]][[j]]
      }
      
      
      # Subset selected features and inner cv splits ------
      X_train_i <- XY_train[split_index_inner[,j], (names(XY_train) %in% c(feature_ids))]  # filter subjects train split + DEA-selected features
      X_val  <- XY_train[-split_index_inner[,j], (names(XY_train) %in% c(feature_ids))] 
      
      Y_train_i <- XY_train[split_index_inner[,j], ] %>% 
        dplyr::select(all_of(target)) %>%
        mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(XY_train[[target]])))))) # target values to 'valid' names: {X0, X1}
      Y_val <- XY_train[-split_index_inner[,j], ] %>% 
        dplyr::select(all_of(target)) %>% 
        mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(XY_train[[target]]))))))
      
      print(paste("Dimensions of inner train set", as.character(dim(X_train_i)[1]), "x", as.character(dim(X_train_i)[2]), sep = " "))
      print(paste("Dimensions of validation set", as.character(dim(X_val)[1]), "x", as.character(dim(X_val)[2]), sep = " "))
      
      ## scale & center data -> should be done within the cv, not before!!
      #preprocess_object <- preProcess(X_train_i, 
      #                                method = c('scale', 'center'))
      #X_train_i <- predict(preprocess_object, X_train_i)
      #X_val <- predict(preprocess_object, X_val)
      
      
      # modelling ------- 
      fitControl <- trainControl(method = "repeatedcv",  # inner cv loop 2
                                 number = 5, # 5-fold cv
                                 repeats = 2,
                                 p = 0.8,
                                 verboseIter = FALSE,
                                 search = "grid",
                                 savePredictions = "final",
                                 classProbs = TRUE,
                                 summaryFunction = twoClassSummary) 
      
      if (model == "glmnet") {
        MODELS[[i]][[model]][[j]] <- caret::train(form = as.formula(paste(target, "~ .")), 
                                                  data = cbind(X_train_i, Y_train_i),  
                                                  method = model,
                                                  trControl = fitControl,
                                                  tuneGrid = expand.grid(alpha = seq(0, 1, length = 15), lambda = parameters),
                                                  metric = "ROC",
                                                  verbose = FALSE,
                                                  preProcess = c('scale', 'center')) # scale & center data
      } else if (model == "deepboost") {
        MODELS[[i]][[model]][[j]] <- caret::train(form = as.formula(paste(target, "~ .")), 
                                                  data = cbind(X_train_i, Y_train_i), 
                                                  method = deepboost_prob, 
                                                  trControl = fitControl,
                                                  metric = "ROC",
                                                  verbose = FALSE,
                                                  preProcess = c('scale', 'center')) # scale & center data
      } else {
        MODELS[[i]][[model]][[j]] <- caret::train(form = as.formula(paste(target, "~ .")), 
                                                  data = cbind(X_train_i, Y_train_i), 
                                                  method = model, 
                                                  trControl = fitControl,
                                                  metric = "ROC",
                                                  verbose = FALSE,
                                                  verbosity = 0, # parameter for xgb-trees models
                                                  preProcess = c('scale', 'center')) # scale & center data
      } 
      
      pred_val = predict(MODELS[[i]][[model]][[j]], newdata = X_val, type = "prob")
      
      # track inner cv performance: mean error, AUC 
      error_pred <- as.factor(predict(MODELS[[i]][[model]][[j]], newdata = X_val))
      error_pred <- mean(ifelse(Y_val %>% pull() != error_pred, 1, 0))
      ERROR[[model]][j, model] <- error_pred
      ERROR[[model]][j, 'fold'] <- j
      
      auc_pred <- pROC::roc(response = as.vector(as.matrix(Y_val)), predictor = pred_val[[2]])$auc
      AUC[[model]][j, model] <- auc_pred
      AUC[[model]][j, 'fold'] <- j
    }
    
    # apply best inner model on outer cv: X_test, Y_test
    best_j <- which.max(AUC[[model]][,1])
    MODELS[[i]][[model]][["best"]] <- MODELS[[i]][[model]][[best_j]]
    pred_test = predict(MODELS[[i]][[model]][[best_j]], newdata = X_test, type = "prob") 
    rownames(pred_test) <- rownames(X_test)
    
    # track outer cv performance: mean error, AUC
    error_pred <- as.factor(predict(MODELS[[i]][[model]][[j]], newdata = X_test))
    print(caret::confusionMatrix(error_pred, Y_test %>% pull(), positive='X1'))
    error_pred <- mean(ifelse(Y_test %>% pull() != error_pred, 1, 0))
    ERROR.OUTER[i, model] <- error_pred
    ERROR.OUTER[i, 'fold'] <- i
    
    auc_pred <- pROC::roc(response = as.vector(as.matrix(Y_test)), predictor = pred_test[[2]], directon = "<")$auc
    AUC.OUTER[i, model] <- auc_pred
    AUC.OUTER[i, 'fold'] <- i
  }
}


# Add {avg, sd} performance of each model in the outer sets (across outer folds)
ERROR.OUTER <- rbind(ERROR.OUTER, c("avg", colMeans(ERROR.OUTER[1:ncol(split_index_outer), 2:ncol(ERROR.OUTER)])))
ERROR.OUTER[, 2:ncol(ERROR.OUTER)] <- lapply(ERROR.OUTER[,2:ncol(ERROR.OUTER)],as.numeric)
ERROR.OUTER <- rbind(ERROR.OUTER, c("sd", colSds(as.matrix(ERROR.OUTER[1:ncol(split_index_outer), 2:ncol(ERROR.OUTER)]))))
ERROR.OUTER[, 2:ncol(ERROR.OUTER)] <- lapply(ERROR.OUTER[,2:ncol(ERROR.OUTER)],as.numeric)
AUC.OUTER <- rbind(AUC.OUTER, c("avg", as.numeric(colMeans(AUC.OUTER[1:ncol(split_index_outer), 2:ncol(AUC.OUTER)]))))
AUC.OUTER[, 2:ncol(AUC.OUTER)] <- lapply(AUC.OUTER[,2:ncol(AUC.OUTER)],as.numeric)
AUC.OUTER <- rbind(AUC.OUTER, c("sd", colSds(as.matrix(AUC.OUTER[1:ncol(split_index_outer), 2:ncol(AUC.OUTER)]))))
AUC.OUTER[, 2:ncol(AUC.OUTER)] <- lapply(AUC.OUTER[,2:ncol(AUC.OUTER)],as.numeric)

print(knitr::kable(ERROR.OUTER, row.names = F))
print(knitr::kable(AUC.OUTER, row.names = F))


# export outputs  --------------------------------------------------------------
readr::write_tsv(ERROR.OUTER, file = OUT.ERROR.FILE) 
readr::write_tsv(AUC.OUTER, file = OUT.AUC.FILE) 
rm(X_train_i, X_val, Y_train_i, Y_val, XY_train, X_test, Y_test, metab.t_train, metab.t_train_i, pheno_train, 
   pheno_train_i, metab.t, pheno, lm_res)





####################
# BUILD FINAL MODEL
####################

# select features for final model: those features repeated at least as many times as outer folds
feature_ids <- as.data.frame(table(unlist(FEAT.SELECTION))) %>% 
  arrange(desc(Freq)) %>% 
  filter(Freq >= ncol(split_index_outer))
print(feature_ids)
feature_ids <- as.character(feature_ids[,1])
metab_4ML <- metab_4ML[, (names(metab_4ML) %in% c(feature_ids, target))] %>%
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(metab_4ML[[target]]))))))

## downsampling not needed

# train final model
rownames(AUC.OUTER) <- NULL
AUC.OUTER <- AUC.OUTER %>%
  column_to_rownames("fold")
model <- names(which.max(AUC.OUTER["avg",]))
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

if (model == "glmnet") {
  f_model <- caret::train(form = as.formula(paste(target, "~ .")), 
                          data = metab_4ML,  
                          method = model,
                          trControl = fitControl,
                          tuneGrid = expand.grid(alpha = seq(0, 1, length = 15), lambda = parameters),
                          metric = "ROC",
                          verbose = FALSE,
                          preProcess = c('scale', 'center')) # scale & center data
} else if (model == "deepboost") {
  f_model <- caret::train(form = as.formula(paste(target, "~ .")), 
                          data = metab_4ML, 
                          method = deepboost_prob, 
                          trControl = fitControl,
                          metric = "ROC",
                          verbose = FALSE,
                          preProcess = c('scale', 'center')) # scale & center data
} else {
  f_model <- caret::train(form = as.formula(paste(target, "~ .")), 
                          data = metab_4ML, 
                          method = model, 
                          trControl = fitControl,
                          metric = "ROC",
                          verbose = FALSE,
                          verbosity = 0, # parameter for xgb-trees models
                          preProcess = c('scale', 'center')) # scale & center data
} 

# final assessment on held-out data
X_hout <- hout_4ML[, !(names(hout_4ML) %in% c(target))][, names(hout_4ML) %in% feature_ids]  
Y_hout <- hout_4ML %>% 
  dplyr::select(all_of(target)) %>% 
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(hout_4ML[[target]]))))))


Y_hat = predict(f_model, newdata = X_hout, type = "prob")
colnames(Y_hat) <- paste(model, colnames(Y_hat), sep = ".")
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
  plot(variable_imp, main=paste(model, "- Variable Importance"))
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
plot(variable_imp_r, main=paste(model, "- Variable Importance (ROC analysis)"))
dev.off()

# output data for shapley values
readr::write_csv(metab_4ML[,-ncol(metab_4ML)] %>% 
                   rownames_to_column(var = "SAMPLE_ID"), file = OUT.SHAPDATA.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

