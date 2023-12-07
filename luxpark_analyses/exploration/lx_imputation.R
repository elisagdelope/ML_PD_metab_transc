# Title: lx_imputation.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs multiple imputation with mice package on batch-normalized metabolomics data
# Usage: Rscript lx_imputation.R
# Data: data from metabolite abundance (batch-normalized metabolomics data) without prior imputation


# GC ----------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(tibble)
library(Hmisc)
library(VIM)
library(mice)
library(missForest)


# I/O --------------------------------------------------------------------------
DIR <- "../data/00-cleansing"
BN.FILE <- file.path(DIR, "batch-norm_data.tsv")
OUT.FILE <- file.path(DIR, "batch-norm_imputed_data_mf.tsv")
OUT_DIR_PLOTS <- DIR


# Data load --------------------------------------------------------------------
data <- vroom(BN.FILE, col_types = cols())
data <- data %>%
  column_to_rownames(var = "SAMPLE_ID")



# Data imputation --------------------------------------------------------------
# 50% cutoff for missing data
pMiss <- function(x){
  sum(is.na(x))/length(x)*100
}

col_missing <- apply(data[,-c(1:2)],2,pMiss) 
row_missing <- apply(data[,-c(1:2)],1,pMiss)
data <- data %>%
  dplyr::select(-any_of(names(col_missing[col_missing > 50]))) %>%
  rownames_to_column(var = "SAMPLE_ID") %>%
  filter(!SAMPLE_ID %in% names(row_missing[row_missing > 50])) %>%
  column_to_rownames(var = "SAMPLE_ID")

# explore missingness
mf <- missForest(data[-c(1:2)], maxiter = 10, verbose = T)
imputed_X <-mf$ximp

imputed_X <- imputed_X %>%
  rownames_to_column(var = "SAMPLE_ID")
data <- data %>%
  rownames_to_column(var = "SAMPLE_ID")

data_complete <- data[,c(1:3)] %>%
  left_join(imputed_X, by = "SAMPLE_ID")

# export output
readr::write_tsv(data_complete, file = OUT.FILE)

# plot histograms
flat_data <- as.vector(as.matrix(data[,-c(1:3)], ncol=1))
flat_data_1 <- flat_data[flat_data < 10]
flat_i_data <- as.vector(as.matrix(imputed_X[,-c(1)], ncol=1))
flat_i_data_1 <- flat_i_data[flat_i_data < 10]
length(flat_data) == length(flat_i_data)

par(mfrow=c(1,2))
plot_name <- "Histograms_distribution_postimputation.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

hist(flat_data_1, breaks = seq(0, 10, 0.05), 
     main = "Distribution without imputed values",
     xlab = "Batch-normalized data")
hist(flat_i_data_1, breaks = seq(0, 10, 0.05), 
     main = "Distribution with imputed values",
     xlab = "Batch-normalized data")

# explore distribution of individual metabolites 
data_6 <- data[, c(4:9)]
data_complete_6 <- data_complete[, c(4:9)]

# plot first 6 columns
par(mfrow=c(1,1))
hist.data.frame(data_6, mtitl = "Distributions of 6 metabolites at random without imputed values")
hist.data.frame(data_complete_6, mtitl = "Distributions of 6 metabolites at random with imputed values")

dev.off()


# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


