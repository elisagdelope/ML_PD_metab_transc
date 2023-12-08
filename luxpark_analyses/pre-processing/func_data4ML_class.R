# Title: func_data4ML_class.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains functions to pre-process data prior to ML modelling

process_data4ML <- function(METAB.FILE, PHENO.FILE, features_varname, OUT_DIR_DATA, e_level, st, target, myseed, export){
  
  # Data load --------------------------------------------------------------------
  # log-transformed peak area & clinical data.
  annotation <- vroom(ANNOTATION.FILE, col_types = cols())
  pheno <- vroom(PHENO.FILE, col_types = cols())
  metab <- vroom(METAB.FILE, col_types = cols()) 
  metab = metab[!duplicated(metab[["SAMPLE_ID"]]), ]
  
  
  
  # Data reformatting for ML -----------------------------------------------------
  if(!all(pheno[[target]] %in% c(0,1))) {
    pheno_4ML <- pheno %>%
      dplyr::select(all_of(c(target, 'SAMPLE_ID'))) %>%
      mutate(!!target := case_when(get(target) == "HC" ~ 0,
                                   get(target) == "PD" ~ 1))
    pheno <- pheno %>% 
      mutate_at(target, factor)
  } else {
    pheno_4ML <- pheno %>%
      dplyr::select(all_of(c(target, 'SAMPLE_ID'))) %>% 
      mutate_at(target, factor) 
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
  
  
  # add target variable information
  metab_4ML <- metab %>%
    dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>%
    inner_join(pheno_4ML, 
               by = "SAMPLE_ID") %>%
    mutate_at(target, factor) %>% 
    column_to_rownames("SAMPLE_ID")

  
  # export pre-processed data 
  if (export == TRUE) {
    if (e_level == "METAB") {
      readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0("data_metab_4ML_", target, ".tsv")))
    } else {
      readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0(e_level, "_", st, "_data_metab_4ML_", target, ".tsv")))
    }
  }
  
  
  
  # create training/held-out set -------------------------------------------------
  set.seed(myseed-1)
  inTraining <- createDataPartition(metab_4ML[[as.character(target)]], p = .85, list = FALSE)
  hout_4ML  <- metab_4ML[-inTraining, ]
  metab_4ML <- metab_4ML[inTraining, ] 
  # export pre-processed data split
  if (export == TRUE) {
    if (e_level == "METAB") {
      readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0("data_cv_metab_4ML_", target, ".tsv")))
      readr::write_tsv(hout_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0("data_test_metab_4ML_", target, ".tsv")))
    } else {
      readr::write_tsv(metab_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0(e_level, "_", st, "_data_cv_metab_4ML_", target, ".tsv")))
      readr::write_tsv(hout_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0(e_level, "_", st, "_data_test_metab_4ML_", target, ".tsv")))
    }
  }
}


process_data4ML_TS <- function(METAB.FILE, pheno_4ML, features_varname, OUT_DIR_DATA, e_level, st, target, temp_feature, myseed, export) {
  
  # Data load --------------------------------------------------------------------
  # log-transformed peak area & clinical data.
  annotation <- vroom(ANNOTATION.FILE, col_types = cols())
  metab <- vroom(METAB.FILE, col_types = cols()) 
  dim_metab <- dim(metab)
  
  metab_4ML <- metab %>%
    inner_join(pheno_4ML, 
               by = "PATIENT_ID") %>%
    mutate_at(target, factor)

  # apply unsupervised filters ---------------------------------------------------
  # remove zero variance features 
  nzv <- nearZeroVar(metab_4ML[, !(names(metab_4ML) %in% c(target))], freqCut = 10, names = TRUE) 
  if (length(nzv) != 0 ) {
    metab_4ML <- metab_4ML %>%
      dplyr::select(-any_of(nzv))
  }
  nzv2 <- nearZeroVar(metab_4ML[, !(names(metab_4ML) %in% c(target))], saveMetrics = TRUE)
  metab_4ML <- metab_4ML %>%
    dplyr::select(-any_of(rownames(nzv2[nzv2$freqRatio > 50,])))
  print(paste("nzv filter:", length(nzv) + length(rownames(nzv2[nzv2$freqRatio > 50,]))))
  
  # remove highly correlated features 
  cor_df = cor(metab_4ML[,-c(1,ncol(metab_4ML))]) # remove patient ID and diagnosis variables
  hc = findCorrelation(cor_df, cutoff=0.85, names = TRUE) 
  hc = sort(hc)
  if (length(hc) > 0) {
    metab_4ML = metab_4ML[,-which(names(metab_4ML) %in% c(hc))]
  }
  print(paste("hc filter:", length(hc)))
  print("NZV and correlation filters successfully applied")
  print(paste((length(hc) + length(nzv) + length(rownames(nzv2[nzv2$freqRatio > 50,]))), "features were removed out of", dim_metab[2]))
  
  metab_4ML <- metab_4ML %>%
    column_to_rownames("PATIENT_ID")
  
  if (e_level == "GENE") {
    readr::write_tsv(expr_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0("data_", temp_feature, "_metab_4ML_", target, ".tsv")))
  } else {
    readr::write_tsv(expr_4ML %>% rownames_to_column(var = "SAMPLE_ID"), file = file.path(OUT_DIR_DATA, paste0(e_level, "_", st, "_", temp_feature, "_data_metab_4ML_", target, ".tsv")))
  }
  
}