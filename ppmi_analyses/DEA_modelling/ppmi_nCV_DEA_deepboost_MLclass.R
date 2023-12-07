# Title: ppmi_nCV_DEA_deepboost_MLclass.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs nested crossvalidation for feature selection (through DEA) + binary deepboost classifier on RNAseq static data (BL).
# Usage (gene level): Rscript ppmi_nCV_DEA_deepboost_MLclass.R 10 -n yes (for 10 features and normalized counts!)
# Usage (aggregated level): Rscript ppmi_nCV_DEA_deepboost_MLclass.R 10 -l gobp -s mean 
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
library(DESeq2)
library(edgeR)
library(argparser, quietly = TRUE)
library(matrixStats)
library(fastshap)
library(deepboost)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-BL-PD"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR_PLOTS <- paste0("../data/", analysis_name , "/01-plots") 
OUT_DIR_PLOTS_SHAP <- paste0(OUT_DIR_PLOTS, "/shap_data") 
PHENO.FILE <- file.path(OUT_DIR, "pheno_BL.tsv")
source("deepboost_prob.R")
target = "DIAGNOSIS" # pass as cmd parameter
myseed = 111

# Add command line arguments
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "k", help = "number of features selected from DEA", default = 10, short = "k")
p <- add_argument(parser = p, arg = "--norm", help = "boolean string (yes/no) whether to use deseq2-flt-normalized counts or not (raw counts)", default = "NO", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
k = as.integer(argv$k) # number of features to select
norm_counts <- toupper(argv$norm)
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level
print(paste("Parameters DEA: k=", k))

if ((e_level == "GENE") & (norm_counts == "NO")) { 
  EXPRS.FILE <- file.path(OUT_DIR, "flt_star_BL.tsv")
  OUT.AUC.FILE <- file.path(OUT_DIR, paste(e_level, paste0(k, "ft"), target, "db_AUC.tsv", sep = "_"))
  OUT.ERROR.FILE <- file.path(OUT_DIR, paste(e_level, paste0(k, "ft"), target, "db_ERROR.tsv", sep = "_"))
  OUT.MODEL.FILE <- file.path(OUT_DIR, paste(e_level, paste0(k, "ft"), target, "db_final_model.rds", sep = "_"))
  OUT.PROBS.FILE <- file.path(OUT_DIR, paste(e_level, paste0(k, "ft"), target, "db_final_probs.tsv", sep = "_"))
  OUT.FPERF.FILE <- file.path(OUT_DIR, paste(e_level, paste0(k, "ft"), target, "db_final_performance.tsv", sep = "_"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, paste0(k, "ft"), target, "db_varImp.pdf", sep = "_"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, paste0(k, "ft"), target, "db_shap_values.csv", sep = "_"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, paste0(k, "ft"), target, "db_shap_data.csv", sep = "_"))
  features_varname <- "GENEID"
  
} else if ((e_level == "GENE") & (norm_counts == "YES")) { # normalized counts at gene level
  EXPRS.FILE <- file.path(OUT_DIR, "flt_norm_star_BL.tsv")
  OUT.AUC.FILE <- file.path(OUT_DIR, paste(e_level, "norm", paste0(k, "ft"), target, "db_AUC.tsv", sep = "_"))
  OUT.ERROR.FILE <- file.path(OUT_DIR, paste(e_level, "norm", paste0(k, "ft"), target, "db_ERROR.tsv", sep = "_"))
  OUT.MODEL.FILE <- file.path(OUT_DIR, paste(e_level, "norm", paste0(k, "ft"), target, "db_final_model.rds", sep = "_"))
  OUT.PROBS.FILE <- file.path(OUT_DIR, paste(e_level, "norm", paste0(k, "ft"), target, "db_final_probs.tsv", sep = "_"))
  OUT.FPERF.FILE <- file.path(OUT_DIR, paste(e_level, "norm", paste0(k, "ft"), target, "db_final_performance.tsv", sep = "_"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, "norm", paste0(k, "ft"), target, "db_varImp.pdf", sep = "_"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", paste0(k, "ft"), target, "db_shap_values.csv", sep = "_"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, "norm", paste0(k, "ft"), target, "db_shap_data.csv", sep = "_"))
  features_varname <- "GENEID"
  
} else { # aggregations
  EXPRS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
  OUT.AUC.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, paste0(k, "ft"), target, "db_AUC.tsv", sep = "_"))
  OUT.ERROR.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, paste0(k, "ft"), target, "db_ERROR.tsv", sep = "_"))
  OUT.MODEL.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, paste0(k, "ft"), target, "db_final_model.rds", sep = "_"))
  OUT.PROBS.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, paste0(k, "ft"), target, "db_final_probs.tsv", sep = "_"))
  OUT.FPERF.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, paste0(k, "ft"), target, "db_final_performance.tsv", sep = "_"))
  OUT.PLOTS.FILE <- file.path(OUT_DIR_PLOTS, paste(e_level, st, paste0(k, "ft"), target, "db_varImp.pdf", sep = "_"))
  OUT.SHAPVALUES.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, paste0(k, "ft"), target, "db_shap_values.csv", sep = "_"))
  OUT.SHAPDATA.FILE <- file.path(OUT_DIR_PLOTS_SHAP, paste(e_level, st, paste0(k, "ft"), target, "db_shap_data.csv", sep = "_"))
  features_varname <- paste0(e_level, "_NAME")
  
  if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier")) | (!file.exists(EXPRS.FILE))) { 
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
rm(pheno_4ML)



# apply unsupervised filters ---------------------------------------------------
# remove zero variance features 
nzv <- nearZeroVar(expr_4ML[, !(names(expr_4ML) %in% c(target))], freqCut = 10 ) 
if (length(nzv) != 0 ) {
  expr_4ML <- expr_4ML[,-nzv]
}

# remove highly correlated features 
cor_df = cor(expr_4ML[,-c(1,ncol(expr_4ML))]) # remove patient ID and diagnosis variables
hc = findCorrelation(cor_df, cutoff=0.85, names = TRUE) 
hc = sort(hc)
expr_4ML = expr_4ML[,-which(names(expr_4ML) %in% c(hc))]
print("NZV and correlation filters successfully applied")
print(paste((length(hc) + length(nzv)), "features were removed out of", nrow(expr)))

# remove genes that didn't pass the filters 
expr <- expr %>% 
  filter(get(features_varname) %in% names(expr_4ML)) %>%  
  column_to_rownames(var = features_varname) 

expr_4ML <- expr_4ML %>% 
  column_to_rownames("PATIENT_ID")



# create training/held-out set -------------------------------------------------
set.seed(myseed-1)
inTraining <- createDataPartition(expr_4ML[[as.character(target)]], p = .85, list = FALSE)
hout_4ML  <- expr_4ML[-inTraining, ]

expr_4ML <- expr_4ML[inTraining, ] 
expr <- expr[, inTraining]
pheno <- pheno[inTraining, ]



# create nested cv -------------------------------------------------------------
set.seed(myseed)
split_index_outer <- createDataPartition(expr_4ML[[target]], p = 0.7, list = FALSE, times = 10)
head(split_index_outer, 10)

# define models & parameters
MODELS = list()
model_list_caret <- c("deepboost") 

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
  expr_train <- expr[, split_index_outer[,i]] 
  pheno_train <- pheno[split_index_outer[,i], ]
  
  XY_train <- expr_4ML[split_index_outer[,i],] # keeps the target var
  X_test  <- expr_4ML[-split_index_outer[,i], !(names(expr_4ML) %in% c(target))]  # without target var
  Y_test <- expr_4ML[-split_index_outer[,i], ] %>% 
    dplyr::select(all_of(target)) %>%  
    mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(XY_train[[target]])))))) # target values to 'valid' names: {X0, X1}
  
  print(paste("Dimensions of outer train set", as.character(dim(XY_train)[1]), "x", as.character(dim(XY_train)[2]-1), sep = " "))
  print(paste("Dimensions of outer test set", as.character(dim(X_test)[1]), "x", as.character(dim(X_test)[2]), sep = " "))
  
  
  # downsampling 
  set.seed(myseed + 1)
  XY_train <- downSample(x=XY_train[,-ncol(XY_train)] %>% rownames_to_column(var = "PATIENT_ID"),
                         y=XY_train[[target]],
                         yname = target) %>%
    arrange(PATIENT_ID) %>% 
    column_to_rownames(var = "PATIENT_ID")
  print("downsampling in inner loop successfully applied:")
  print(table(XY_train[[target]]))
  
  # inner loop partitions
  set.seed(myseed + 2)
  split_index_inner <- createDataPartition(XY_train[[target]], p = 0.7, list = FALSE, times = 5)
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
        
        if ((st == "pathifier") & (e_level != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
          expr_train_i <- round(expr_train[, split_index_inner[,j]] * 100, digits = 0) # *100 AND remove decimal digits -> not allowed for the DEA 
        } else {
          expr_train_i <- round(expr_train[, split_index_inner[,j]], digits = 0) # remove decimal digits (from the normalization) -> not allowed for the DEA 
        }
        expr_train_i.deseq <- DESeq2::DESeqDataSetFromMatrix(expr_train_i, 
                                                             colData = pheno_train_i, 
                                                             design = ~ AGE + GENDER + DIAGNOSIS)
        expr_train_i.dgelist <- DGEList(counts(expr_train_i.deseq), group = colData(expr_train_i.deseq)$DIAGNOSIS) # The function accesses the group factor to compute the minimum group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced
        gene.flt <- edgeR::filterByExpr(expr_train_i.dgelist, min.count = 10)
        expr_train_i.flt <- expr_train_i.deseq[gene.flt, ]
        
        if ((e_level != "GENE") | ((e_level == "GENE") & (norm_counts == "YES"))) {
          # deseq ATTENTION: DATA ALREADY NORMALIZED
          sizeFactors_nonnorm <- rep(1, ncol(expr_train_i.flt)) # create a vector of 1s to avoid normalization
          names(sizeFactors_nonnorm) <- colnames(expr_train_i.flt)
          sizeFactors(expr_train_i.flt) <- sizeFactors_nonnorm
          dds <- estimateDispersions(expr_train_i.flt)
          dds <- nbinomWaldTest(dds, maxit=500)
          res <- results(dds, contrast = c("DIAGNOSIS", "PD", "HC"))
        } else {
          dds <- DESeq(expr_train_i.flt) # complete DESEQ process (for not-notmalized data)
          res <- results(dds, contrast = c("DIAGNOSIS", "PD", "HC"))
        }
        
        res <- as.data.frame(res) %>% 
          rownames_to_column(var = features_varname) %>%
          filter(!is.na(padj)) %>%
          arrange(pvalue) %>%
          slice_head(n = k)
        
        feature_ids <- res[[1]] 
        FEAT.SELECTION[[i]][[j]] <- res[[1]]
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
      
      ## Still need to preprocess each new set! TO BE REMOVED
      #preprocess_object <- preProcess(X_train_i, 
      #                                method = c('scale', 'center'))
      #X_train_i <- predict(preprocess_object, X_train_i)
      #X_val <- predict(preprocess_object, X_val)
      
      
      # modelling ------- 
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
    
      MODELS[[i]][[model]][[j]] <- caret::train(form = as.formula(paste(target, "~ .")), 
                                                  data = cbind(X_train_i, Y_train_i), 
                                                  method = deepboost_prob, 
                                                  trControl = fitControl,
                                                  metric = "ROC",
                                                  verbose = FALSE) # parameter for xgb-trees models
      
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
    print(error_pred)
    ERROR.OUTER[i, model] <- error_pred
    ERROR.OUTER[i, 'fold'] <- i
    
    auc_pred <- pROC::roc(response = as.vector(as.matrix(Y_test)), predictor = pred_test[[2]], directon = "<")$auc
    print(auc_pred)
    AUC.OUTER[i, model] <- auc_pred
    AUC.OUTER[i, 'fold'] <- i
  }
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
rm(X_train_i, X_val, Y_train_i, Y_val, XY_train, X_test, Y_test, expr_train, expr_train_i, pheno_train, 
   pheno_train_i, expr, pheno, res, expr_train_i.deseq, expr_train_i.dgelist, expr_train_i.flt)



####################
# BUILD FINAL MODEL
####################

# select features for final model: those features repeated at least as many times as outer folds
feature_ids <- as.data.frame(table(unlist(FEAT.SELECTION))) %>% 
  arrange(desc(Freq)) %>% 
  filter(Freq >= ncol(split_index_outer))
print(feature_ids)
feature_ids <- as.character(feature_ids[,1])
expr_4ML <- expr_4ML[, (names(expr_4ML) %in% c(feature_ids, target))] %>%
  mutate(across(all_of(target), ~factor(.x, labels = make.names(sort(unique(expr_4ML[[target]]))))))

# downsampling 
set.seed(myseed + 3)
expr_4ML <- downSample(x=expr_4ML[,-ncol(expr_4ML)] %>% rownames_to_column(var = "PATIENT_ID"),
                       y=expr_4ML[[target]],
                       yname = target) %>%
  arrange(PATIENT_ID) %>% 
  column_to_rownames(var = "PATIENT_ID")
print("downsampling in final training successfully applied:")
print(table(expr_4ML[[target]]))

# train final model
mymodel <- "deepboost"
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
                          method = deepboost_prob, 
                          trControl = fitControl,
                          metric = "ROC",
                          verbose = FALSE) # parameter for xgb-trees models


# final assessment on held-out data
X_hout <- hout_4ML[, !(names(hout_4ML) %in% c(target))][, names(hout_4ML) %in% feature_ids]  
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
plot(variable_imp, main=paste(mymodel, "- Variable Importance"))
variable_imp_r <- caret::varImp(f_model, useModel = FALSE, scale = FALSE) #from caret 
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
