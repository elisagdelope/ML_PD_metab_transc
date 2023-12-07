# Title: lx_plot_distributions.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script explores and generates distribution plots for the full data matrix of metabolomics data + distribution of 6 individual metabolites at random
# Usage: R lx_plot_distributions.R
# Data: data from metabolite abundance pre and post imputation


# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)


# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(tibble)
library(Hmisc)



# I/O --------------------------------------------------------------------------
BN.FILE <- "../data/00-cleansing/batch-norm_data.tsv"
BN_IMPUTED.FILE <- "../data/00-cleansing/batch-norm_imputed_data.tsv"
OUT_DIR_PLOTS <- "../data/00-cleansing/"


# Data load --------------------------------------------------------------------
data <- vroom(BN.FILE, col_types = cols())
i_data <- vroom(BN_IMPUTED.FILE, col_types = cols())



# exploration ------------------------------------------------------------------
# full df
flat_data <- as.vector(as.matrix(data[,-c(1:3)], ncol=1))
flat_data_1 <- flat_data[flat_data < 10]
flat_i_data <- as.vector(as.matrix(i_data[,-c(1:3)], ncol=1))
flat_i_data_1 <- flat_i_data[flat_i_data < 10]

# plot
par(mfrow=c(1,2))
plot_name <- "Histograms_distribution.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

hist(flat_data_1, breaks = seq(0, 10, 0.05), 
     main = "Distribution of full data matrix without imputed values",
     xlab = "Batch-normalized data")
hist(flat_i_data_1, breaks = seq(0, 10, 0.05), 
     main = "Distribution of full data matrix with imputed values",
     xlab = "Batch-normalized data")

# explore distribution of individual metabolites 
data_6 <- data[, c(4:9)]
i_data_6 <- i_data[, c(4:9)]

# plot first 6 columns
par(mfrow=c(1,1))
hist.data.frame(data_6, mtitl = "Distributions of 6 metabolites at random without imputed values")
hist.data.frame(i_data_6, mtitl = "Distributions of 6 metabolites at random with imputed values")

dev.off()



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


