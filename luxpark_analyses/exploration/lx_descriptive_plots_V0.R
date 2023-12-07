# Title: lx_descriptive_plots_V0.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates descriptive stats plots (PCoA based on different distance metrics & PCA, dendrograms, distribution per visits, trajectories,...)
# Usage: R lx_descriptive_plots_V0.R
# Data: data from metabolite abundance (filtered + transformed) + pheno (filtered for training)

# GC ----------------------------------------------------------------------
rm(list = ls())
gc(T)


# Packages ----------------------------------------------------------------
#library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(FactoMineR)
library(ggplot2)
library(tibble)
library(MASS)
library(RColorBrewer)
library(viridis)
library(sparcl)


# I/O ---------------------------------------------------------------------
analysis_name <- "01-dea-V0-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
EXPRESSION.FILE <- paste0("../data/", analysis_name, "/02-outfiles/log_transformed_V0.tsv")
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv")


# Main --------------------------------------------------------------------
if (!dir.exists(OUT_DIR_PLOTS)) {
  dir.create(OUT_DIR_PLOTS, recursive = T)
}


# Data load --------------------------------------------------------------
expression <- vroom(EXPRESSION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 


# PCA -----------------------------------------------------------
coul=rev(rainbow(6))

# remove low variance features
nzv <- nearZeroVar(expression[,-c(1:3)], names = TRUE)
expression <- expression %>%
  dplyr::select(-any_of(nzv))

data <- expression[, -c(1:2)] %>%
  column_to_rownames("SAMPLE_ID")


# PCoA data euclidean distance (patients)
dist_matrix <- dist(data, method = "euclidean")
pcoa <- stats::cmdscale(dist_matrix)
pcoa_pheno <- pcoa %>%
  as.data.frame() %>%
  rownames_to_column(var = 'SAMPLE_ID') %>%
  left_join(pheno %>% dplyr::select(PATIENT_ID, VISIT, SAMPLE_ID, DIAGNOSIS, GENDER), by= "SAMPLE_ID")

# PCoA data spearman distance (patients)
dist_matrix_spearman <- as.dist(1 - abs(cor(t(data), method = "spearman")))
pcoa_spearman <- stats::cmdscale(dist_matrix_spearman)
pcoa_spearman_pheno <- pcoa_spearman %>%
  as.data.frame() %>%
  rownames_to_column(var = 'SAMPLE_ID') %>%
  left_join(pheno %>% dplyr::select(PATIENT_ID, VISIT, SAMPLE_ID, DIAGNOSIS, GENDER), by= "SAMPLE_ID")

# PCoA data pearson distance (patients)
dist_matrix_pearson <- as.dist(1 - abs(cor(t(data), method = "pearson")))
pcoa_pearson <- stats::cmdscale(dist_matrix_pearson)
pcoa_pearson_pheno <- pcoa_pearson %>%
  as.data.frame() %>%
  rownames_to_column(var = 'SAMPLE_ID') %>%
  left_join(pheno %>% dplyr::select(PATIENT_ID, VISIT, SAMPLE_ID, DIAGNOSIS, GENDER), by= "SAMPLE_ID")

# PCA (prcomp function)
pca <- prcomp(data, center = TRUE, scale. = TRUE)
pca_pheno <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column(var = 'SAMPLE_ID') %>%
  left_join(pheno %>% dplyr::select(PATIENT_ID, VISIT, SAMPLE_ID, DIAGNOSIS, GENDER), by= "SAMPLE_ID")


data_pheno <- merge(data, pheno %>% dplyr::select(PATIENT_ID, VISIT, SAMPLE_ID, DIAGNOSIS, GENDER), by.x = "row.names", by.y ="SAMPLE_ID") %>%
  tibble::column_to_rownames(var = "Row.names") %>%
  mutate(DIAGNOSIS_NUM = case_when(DIAGNOSIS == "HC" ~ 1,
                                   DIAGNOSIS == "PD" ~ 2))

# PCA (FactoMineR)
pca2 <- FactoMineR::PCA(data_pheno[,1:(ncol(data_pheno) - 5)])


# PLOT PCoA, PCA 
plot_name <- "Descriptive_statistics_PCA_PCoA.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

pcoa_pheno %>%
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = GENDER)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (euclidean distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_spearman_pheno %>%
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = GENDER)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (spearman corr distance) of all samples")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_pearson_pheno %>%
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = GENDER)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (pearson corr distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pca_pheno %>%
  ggplot(aes(PC1, PC2, colour =  DIAGNOSIS, shape = GENDER)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCA of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

# additional pca plot
plot(as.numeric(pca2$ind$coord[,1]),as.numeric(pca2$ind$coord[,2]),
     col = coul[data_pheno$DIAGNOSIS_NUM], 
     # ylim=c(-60,60),xlim=c(-50,50),
     pch = c(15, 16, 17, 18)[as.numeric(factor(data_pheno$VISIT))],
     main = "PCA of all samples", 
     xlab = paste("Dim1: ",round(pca2$eig[1,2],2),"%",sep = ""),
     ylab = paste("Dim2: ",round(pca2$eig[2,2],2),"%",sep = ""))
abline(v = 0,
       h = 0,
       lty = 3)
legend("topleft",
       legend = c(levels(factor(data_pheno$DIAGNOSIS)),levels(factor(data_pheno$VISIT))),
       pch = c(rep(16, length(levels(factor(data_pheno$DIAGNOSIS)))), c(15, 16, 17, 18)), 
       col = c(levels(factor(coul[data_pheno$DIAGNOSIS_NUM])), rep("#000000", length(levels(factor(data_pheno$VISIT))))), 
       bty = 'o', bg = 'lightgrey', x.intersp = 0.75, title = "Diagnosis / timepoints", horiz = TRUE, text.width = 10)
dev.off()




# CLUSTERING ----------------------------------------------------------------------

plot_name <- "Descriptive_statistics_dendrograms.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

distance2=dist(data,method="euclidean")
hier2=hclust(distance2,method="ward.D2") 
ColorDendrogram(hier2,y=coul[data_pheno$DIAGNOSIS_NUM],labels=rownames(data_pheno),
                main="Dendrogram (Ward's method) of metabolite abundance by diagnosis",xlab="sample",sub="",branchlength = 300000)
dev.off()



