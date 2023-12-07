# Title: lx_descriptive_plots.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates descriptive stats plots (PCoA based on different distance metrics & PCA, dendrograms, distribution per visits, trajectories,...)
# Usage: R lx_descriptive_plots.R
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
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
#OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
EXPRESSION.FILE <- "../data/00-cleansing/log_transformed_data.tsv"
PHENO.FILE <- "../data/00-cleansing/lx_pheno.tsv"


# Main --------------------------------------------------------------------
if (!dir.exists(OUT_DIR_PLOTS)) {
  dir.create(OUT_DIR_PLOTS, recursive = T)
}


# Data load --------------------------------------------------------------
expression <- vroom(EXPRESSION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = cols())



# Filter out metabolites with low variance -------------------------------------

#var_filter = function(X, filtsize=1000)
#{
#  filtsize <- min(ncol(X), as.numeric(filtsize))
#  variances <- apply(X, 2, var) # apply variance on columns
#  o <- order(variances, decreasing = TRUE)
#  Xfilt <- (X[,o])[,1:filtsize]
#  return(Xfilt)
#}
# filter the expression matrices to only retain the top 1000 metabolites with the highest variance
#expression = as.data.frame(cbind(expression[,c(1:3)], var_filter(expression[,-c(1:3)], filtsize = 1000)))
nzv <- nearZeroVar(expression[,-c(1:3)], names = TRUE)
expression <- expression %>%
  dplyr::select(-any_of(nzv))

# remove visits with <10 patients:
v <- names(which(table(expression$VISIT) < 20))
expression <- expression %>%
  filter(!VISIT %in% v)
  

# PCA -----------------------------------------------------------
coul=rev(rainbow(6))
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
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = VISIT, size= GENDER)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (euclidean distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_spearman_pheno %>%
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = VISIT, size= GENDER)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (spearman corr distance) of all samples")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_pearson_pheno %>%
  ggplot(aes(V1, V2, colour = DIAGNOSIS, shape = VISIT, size= GENDER)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (pearson corr distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pca_pheno %>%
  ggplot(aes(PC1, PC2, colour =  DIAGNOSIS, shape = VISIT, size= GENDER)) +
  geom_point() + 
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
ColorDendrogram(hier2,y=coul[as.numeric(factor(data_pheno$VISIT))],labels=as.character(data_pheno$DIAGNOSIS),
                main="Dendrogram (Ward's method) of metabolite abundance by timepoint and diagnosis",xlab="sample",sub="",branchlength = 300000)
dev.off()



# PLOT DISTRIBUTIONS -----------------------------------------------------------

# BOXPLOTS
plot_name <- "Boxplots_genes.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)
# plot boxpolots for 100 genes 
for(i in 1:100){
  boxplot(data_pheno[,i]~data_pheno$VISIT, main=colnames(data_pheno)[i], outline = F, ylim=c(min(data_pheno[,i]),max(data_pheno[,i])),border = coul)
  points(jitter(rep(1,nrow(data_pheno[data_pheno$VISIT=="V0",])),3),data_pheno[data_pheno$VISIT=="V0",][,i],col=coul[1],pch=16)
  points(jitter(rep(2,nrow(data_pheno[data_pheno$VISIT=="V1",])),3),data_pheno[data_pheno$VISIT=="V1",][,i],col=coul[2],pch=16)
  points(jitter(rep(3,nrow(data_pheno[data_pheno$VISIT=="V2",])),3),data_pheno[data_pheno$VISIT=="V2",][,i],col=coul[3],pch=16)
  points(jitter(rep(4,nrow(data_pheno[data_pheno$VISIT=="V3",])),3),data_pheno[data_pheno$VISIT=="V3",][,i],col=coul[4],pch=16)
}
dev.off()


plot_name <- "Violins_visits.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

exprs_pheno <- expression %>%
  pivot_longer(cols = stringr::str_subset(colnames(expression), "M[0-9]"), names_to = "METAB", values_to = "COUNTS") %>%
  left_join(pheno, by = c("SAMPLE_ID", "PATIENT_ID", "VISIT"))

#violin + points
exprs_pheno %>%
  ggplot(aes(x = VISIT,
             y = COUNTS)) +
  geom_violin(trim = TRUE) +
  geom_point() +
  labs(title="Distribution of metabolite abundance by timepoints",
       x ="Time points", 
       y = "Metab abundance")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

# violin + boxplot
exprs_pheno %>%
  ggplot(aes(x = VISIT,
             y = COUNTS,
             fill = DIAGNOSIS)) +
  geom_violin(position = "dodge", trim = FALSE) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  labs(title = "Distribution of metabolite abundance by timepoints",
       x = "Time points", 
       y = "Metab abundance")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))


# violin + mean
exprs_pheno %>%
  ggplot(aes(x = VISIT,
             y = COUNTS,
             fill = DIAGNOSIS)) +
  geom_violin(position = "dodge", trim = FALSE)  +
  stat_summary(fun=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(title="Distribution of metabolite abundance and mean by diagnosis and timepoints",
       x ="Time points", 
       y = "Metab abundance") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))


# violin + boxplot + sample size (number of patients per visit)
sample_size = exprs_pheno %>% 
  group_by(VISIT) %>% 
  summarize(num=n()/length(unique(exprs_pheno$METAB)))
exprs_pheno %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(VISIT, "\n", "n=", num)) %>% 
  ggplot(aes(x = myaxis,
             y = COUNTS,
             fill = VISIT)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete=TRUE) +
  theme_linedraw() +
  labs(title="Distribution of metabolite abundance and sample size by diagnosis and timepoints",
       x ="Time points", 
       y = "Metab abundance") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

dev.off()

# PLOT TRAJECTORIES -----------------------------------------------
palette <- brewer.pal(6, "Set1") 

genes_traj <- exprs_pheno %>% 
  group_by(METAB, VISIT) %>% 
  summarise(median_patients = median(COUNTS)) %>%
  pivot_wider(names_from=VISIT, values_from=median_patients) %>%
  column_to_rownames(var = "METAB") %>% 
  drop_na() # remove genes that have NAs at one or more timepoints

patients_traj <- exprs_pheno %>% 
  group_by(PATIENT_ID, VISIT) %>% 
  summarise(median_patients = median(COUNTS)) %>%
  pivot_wider(names_from=VISIT, values_from=median_patients) %>%
  column_to_rownames(var = "PATIENT_ID") %>% 
  drop_na() # remove patients that have NAs at one or more timepoints

# patient-phenotype dataset for color code:
pat_pheno <- exprs_pheno %>%
  group_by(PATIENT_ID) %>%
  summarise(DIAGNOSIS = unique(DIAGNOSIS),
            GENDER = unique(GENDER))



## output genes and patient trajectories
#readr::write_tsv(genes_traj %>% 
#                   rownames_to_column(var = "geneid"), file = file.path(OUT_DIR, "patients_trajectory.tsv"))
#
#readr::write_tsv(patients_traj %>% 
#                   rownames_to_column(var = "patientid"), file = file.path(OUT_DIR, "genes_trajectory.tsv"))


plot_name <- "TS_trajectories_visits.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

# time series plot(for each gene: median exprs per visit)
# plot only 100-200 genes, otherwise plot is overwhelming
parcoord(genes_traj[1:200,], 
         col=palette,
         main = "Median trajectories of counts for 200 genes")

# time series plot(for each patient: median exprs per visit)
# can't color by diagnosis since all are PD
parcoord(patients_traj, 
         main = "Median trajectories of gene counts for 200 random patients by diagnosis")
legend("topleft",legend=levels(factor(pat_pheno$DIAGNOSIS)), pch = 16, bty = 'n', pt.bg = 'lightgrey', x.intersp = 0.5, horiz = FALSE, title = "Diagnosis")

# color by Gender
my_colors <- palette[as.numeric(factor(pat_pheno$GENDER))]
parcoord(patients_traj[1:200,], 
         col=my_colors,
         main = "Median trajectories of gene counts for 200 random patients by gender")
legend("topleft",legend=levels(factor(pat_pheno$GENDER)), pch = 16, col=my_colors, bty = 'n', pt.bg='lightgrey', x.intersp = 0.5, horiz = FALSE, title = "Gender")



# time series plot(for diagnosis: median exprs of all genes on all patients per visit)
diagnosis_traj <- exprs_pheno %>% 
  group_by(DIAGNOSIS, VISIT) %>% 
  summarise(median_diagnosis = median(COUNTS)) 

my_colors <- palette[as.numeric(factor(rownames(diagnosis_traj)))]
diagnosis_traj %>%
  ggplot(aes(x = VISIT,
             y = median_diagnosis,
             group = DIAGNOSIS,
             color = DIAGNOSIS)) +
  geom_line() +
  labs(title="Median trajectories of gene counts on all patients by diagnosis",
       x ="Time points", 
       y = "Median gene counts") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

dev.off()


