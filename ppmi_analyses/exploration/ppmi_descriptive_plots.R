# Title: ppmi_descriptive_plots.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates descriptive stats plots (PCoA based on different distance metrics & PCA, dendrograms, distribution per visits, trajectories,...)
# Usage: R ppmi_descriptive_plots.R
# Data: data from gene expression (filtered + transformed) + pheno (filtered for training)

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
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
EXPRESSION.FILE <- file.path(OUT_DIR, "flt_transf_star_all_training_TS.tsv") 
PHENO.FILE <- file.path(OUT_DIR, "ppmi_pheno_dsq_training.tsv")


# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) || (!dir.exists(OUT_DIR_PLOTS))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
}


# Data load --------------------------------------------------------------
expression <- vroom(EXPRESSION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = cols())



# Filter out genes with low variance -------------------------------------

var_filter = function(X, filtsize=1000)
{
  filtsize <- min(nrow(X), as.numeric(filtsize))
  variances <- apply(X, 1, var) # apply variance on rows
  o <- order(variances, decreasing = TRUE)
  Xfilt <- (X[o,])[1:filtsize,]
  return(Xfilt)
}


# filter the expression matrices to only retain the top 2000 genes with the highest variance
expression = var_filter(expression %>%
                         column_to_rownames(var = 'geneid'), filtsize = 2000)


# PCA -----------------------------------------------------------
coul=rev(rainbow(6))
data <- t(expression) 

# PCoA data euclidean distance (patients)
dist_matrix <- dist(data, method = "euclidean")
pcoa <- stats::cmdscale(dist_matrix)
pcoa_pheno <- pcoa %>%
              as.data.frame() %>%
              rownames_to_column(var = 'geneid') %>%
              dplyr::mutate(visit = str_extract(rownames(pcoa), "[A-Z]+\\d*$"))  %>%
              left_join(pheno %>% dplyr::select(patient_visit_id, Diagnosis, Gender), by= c("geneid" = "patient_visit_id"))

# PCoA data spearman distance (patients)
dist_matrix_spearman <- as.dist(1 - abs(cor(t(data), method = "spearman")))
pcoa_spearman <- stats::cmdscale(dist_matrix_spearman)
pcoa_spearman_pheno <- pcoa_spearman %>%
  as.data.frame() %>%
  rownames_to_column(var = 'geneid') %>%
  dplyr::mutate(visit = str_extract(rownames(pcoa_spearman), "[A-Z]+\\d*$"))  %>%
  left_join(pheno %>% dplyr::select(patient_visit_id, Diagnosis, Gender), by= c("geneid" = "patient_visit_id"))

# PCoA data pearson distance (patients)
dist_matrix_pearson <- as.dist(1 - abs(cor(t(data), method = "pearson")))
pcoa_pearson <- stats::cmdscale(dist_matrix_pearson)
pcoa_pearson_pheno <- pcoa_pearson %>%
  as.data.frame() %>%
  rownames_to_column(var = 'geneid') %>%
  dplyr::mutate(visit = str_extract(rownames(pcoa_pearson), "[A-Z]+\\d*$"))  %>%
  left_join(pheno %>% dplyr::select(patient_visit_id, Diagnosis, Gender), by= c("geneid" = "patient_visit_id"))

# PCA (prcomp function)
pca <- prcomp(data, center = TRUE, scale. = TRUE)
pca_pheno <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column(var = 'geneid') %>%
  dplyr::mutate(visit = str_extract(rownames(pcoa), "[A-Z]+\\d*$"))  %>%
  left_join(pheno %>% dplyr::select(patient_visit_id, Diagnosis, Gender), by= c("geneid" = "patient_visit_id"))

data_pheno <- merge(data, pheno %>% dplyr::select(patient_visit_id, visit, Diagnosis, Gender), by.x = "row.names", by.y ="patient_visit_id") %>%
  tibble::column_to_rownames(var = "Row.names") %>%
  mutate(Diagnosis_num = case_when(Diagnosis == "HC" ~ 1,
                                   Diagnosis == "PD" ~ 2))

# PCA (FactoMineR)
pca2 <- FactoMineR::PCA(data_pheno[,1:(ncol(data_pheno) - 4)])



# PLOT PCoA, PCA 
plot_name <- "Descriptive_statistics_PCA_PCoA.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

pcoa_pheno %>%
  ggplot(aes(V1, V2, colour = Diagnosis, shape = visit, size= Gender)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (euclidean distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_spearman_pheno %>%
  ggplot(aes(V1, V2, colour = Diagnosis, shape = visit, size= Gender)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (spearman corr distance) of all samples")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

pcoa_pearson_pheno %>%
  ggplot(aes(V1, V2, colour = Diagnosis, shape = visit, size= Gender)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCoA (pearson corr distance) of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))


pca_pheno %>%
  ggplot(aes(PC1, PC2, colour = Diagnosis, shape = visit, size= Gender)) +
  geom_point() + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(title = "PCA of all samples") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))
  

# additional pca plot
plot(as.numeric(pca2$ind$coord[,1]),as.numeric(pca2$ind$coord[,2]),
     col = coul[data_pheno$Diagnosis_num], 
    # ylim=c(-60,60),xlim=c(-50,50),
     pch = c(15, 16, 17, 18)[as.numeric(factor(data_pheno$visit))],
     main = "PCA of all samples", 
     xlab = paste("Dim1: ",round(pca2$eig[1,2],2),"%",sep = ""),
     ylab = paste("Dim2: ",round(pca2$eig[2,2],2),"%",sep = ""))
abline(v = 0,
       h = 0,
       lty = 3)
legend("topleft",
       legend = c(levels(factor(data_pheno$Diagnosis)),levels(factor(data_pheno$visit))),
       pch = c(rep(16, length(levels(factor(data_pheno$Diagnosis)))), c(15, 16, 17, 18)), 
       col = c(levels(factor(coul[data_pheno$Diagnosis_num])), rep("#000000", length(levels(factor(data_pheno$visit))))), 
       bty = 'o', bg = 'lightgrey', x.intersp = 0.75, title = "Diagnosis / timepoints", horiz = TRUE, text.width = 10)
dev.off()


# X <- matrix(rnorm(1000), nrow=10)
# #d <- dist(X)
# #pcoa.res <- cmdscale(d)
# pca.res <- prcomp(X, center = TRUE, scale. = TRUE)
# pca2.res <- FactoMineR::PCA(X)
# plot(c(pca2.res$ind$coord[,1]), c(pca.res$x[, 1]))
# plot(c(pca2.res$ind$coord[,2]), c(pca.res$x[, 2]))
# 
# plot(c(pca2$ind$coord[,1]), c(pca$x[, 1]))
# plot(c(pca2$ind$coord[,2]), c(pca$x[, 2]))




# CLUSTERING ----------------------------------------------------------------------

plot_name <- "Descriptive_statistics_dendrograms.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

distance2=dist(data,method="euclidean")
hier2=hclust(distance2,method="ward.D2") 
ColorDendrogram(hier2,y=coul[data_pheno$Diagnosis_num],labels=rownames(data_pheno),
                main="Dendrogram (Ward's method) of gene counts by diagnosis",xlab="sample",sub="",branchlength = 300000)
ColorDendrogram(hier2,y=coul[as.numeric(factor(data_pheno$visit))],labels=as.character(data_pheno$Diagnosis),
                main="Dendrogram (Ward's method) of gene counts by timepoint and diagnosis",xlab="sample",sub="",branchlength = 300000)
dev.off()



# PLOT DISTRIBUTIONS -----------------------------------------------------------

# BOXPLOTS
plot_name <- "Boxplots_genes.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)
# plot boxpolots for 100 genes 
for(i in 1:100){
  boxplot(data_pheno[,i]~data_pheno$visit, main=colnames(data_pheno)[i], outline = F, ylim=c(min(data_pheno[,i]),max(data_pheno[,i])),border = coul)
  points(jitter(rep(1,nrow(data_pheno[data_pheno$visit=="BL",])),3),data_pheno[data_pheno$visit=="BL",][,i],col=coul[1],pch=16)
  points(jitter(rep(2,nrow(data_pheno[data_pheno$visit=="V04",])),3),data_pheno[data_pheno$visit=="V04",][,i],col=coul[2],pch=16)
  points(jitter(rep(3,nrow(data_pheno[data_pheno$visit=="V06",])),3),data_pheno[data_pheno$visit=="V06",][,i],col=coul[3],pch=16)
  points(jitter(rep(4,nrow(data_pheno[data_pheno$visit=="V08",])),3),data_pheno[data_pheno$visit=="V08",][,i],col=coul[4],pch=16)
}
dev.off()



plot_name <- "Violins_visits.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

exprs_pheno <- expression %>%
  rownames_to_column(var = "geneid") %>%
  pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "counts") %>%
  left_join(pheno, by = c("sample" = "patient_visit_id"))

#violin + points
exprs_pheno %>%
  ggplot(aes(x = visit,
          y = counts)) +
  geom_violin(trim = TRUE) +
  geom_point() +
  labs(title="Distribution of gene counts by timepoints",
       x ="Time points", 
       y = "Gene counts")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

# violin + boxplot
exprs_pheno %>%
  ggplot(aes(x = visit,
             y = counts,
             fill = Diagnosis)) +
  geom_violin(position = "dodge", trim = FALSE) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  labs(title = "Distribution of gene counts by diagnosis and timepoints",
       x = "Time points", 
       y = "Gene counts")  + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

  
# violin + mean
exprs_pheno %>%
  ggplot(aes(x = visit,
             y = counts,
             fill = Diagnosis)) +
  geom_violin(position = "dodge", trim = FALSE)  +
  stat_summary(fun=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(title="Distribution of gene counts and mean by diagnosis and timepoints",
       x ="Time points", 
       y = "Gene counts") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))


# violin + boxplot + sample size (number of patients per visit)
sample_size = exprs_pheno %>% 
              group_by(visit) %>% 
              summarize(num=n()/length(unique(exprs_pheno$geneid)))
exprs_pheno %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(visit, "\n", "n=", num)) %>% 
  ggplot(aes(x = myaxis,
             y = counts,
             fill = visit)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete=TRUE) +
  theme_linedraw() +
  labs(title="Distribution of gene counts and sample size by diagnosis and timepoints",
       x ="Time points", 
       y = "Gene counts") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

dev.off()

# PLOT TRAJECTORIES -----------------------------------------------
palette <- brewer.pal(6, "Set1") 

genes_traj <- exprs_pheno %>% 
  group_by(geneid, visit) %>% 
  summarise(median_patients = median(counts)) %>%
  pivot_wider(names_from=visit, values_from=median_patients) %>%
  column_to_rownames(var = "geneid") %>% 
  drop_na() # remove genes that have NAs at one or more timepoints

patients_traj <- exprs_pheno %>% 
  group_by(GSEID, visit) %>% 
  summarise(median_patients = median(counts)) %>%
  pivot_wider(names_from=visit, values_from=median_patients) %>%
  column_to_rownames(var = "GSEID") %>% 
  drop_na() # remove patients that have NAs at one or more timepoints

# patient-phenotype dataset for color code:
pat_pheno<- exprs_pheno %>%
  group_by(GSEID) %>%
  summarise(Diagnosis = unique(Diagnosis),
            Gender = unique(Gender))



# output genes and patient trajectories
readr::write_tsv(genes_traj %>% 
                   rownames_to_column(var = "geneid"), file = file.path(OUT_DIR, "patients_trajectory.tsv"))

readr::write_tsv(patients_traj %>% 
                   rownames_to_column(var = "patientid"), file = file.path(OUT_DIR, "genes_trajectory.tsv"))


plot_name <- "TS_trajectories_visits.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)

# time series plot(for each gene: median exprs per visit)
# plot only 100-200 genes, otherwise plot is overwhelming
parcoord(genes_traj[1:200,], 
         col=palette,
         main = "Median trajectories of counts for 200 genes")

# time series plot(for each patient: median exprs per visit)
# color by diagnosis
ispd <- ifelse(pat_pheno$Diagnosis =="PD", "red", "blue")
parcoord(patients_traj[1:200,], 
         col=ispd,
         main = "Median trajectories of gene counts for 200 random patients by diagnosis")
legend("topleft",legend=levels(factor(pat_pheno$Diagnosis)), pch = 16, col = ispd, bty = 'n', pt.bg = 'lightgrey', x.intersp = 0.5, horiz = FALSE, title = "Diagnosis")

# color by Gender
my_colors <- palette[as.numeric(factor(pat_pheno$Gender))]
parcoord(patients_traj[1:200,], 
         col=my_colors,
         main = "Median trajectories of gene counts for 200 random patients by gender")
legend("topleft",legend=levels(factor(pat_pheno$Gender)), pch = 16, col=my_colors, bty = 'n', pt.bg='lightgrey', x.intersp = 0.5, horiz = FALSE, title = "Gender")



# time series plot(for diagnosis: median exprs of all genes on all patients per visit)
diagnosis_traj <- exprs_pheno %>% 
  group_by(Diagnosis, visit) %>% 
  summarise(median_diagnosis = median(counts)) 

my_colors <- palette[as.numeric(factor(rownames(diagnosis_traj)))]
diagnosis_traj %>%
  ggplot(aes(x = visit,
             y = median_diagnosis,
             group = Diagnosis,
             color = Diagnosis)) +
  geom_line() +
  labs(title="Median trajectories of gene counts on all patients by diagnosis",
        x ="Time points", 
        y = "Median gene counts") + 
  theme(plot.title = element_text(color = "black", size = 12, face = "bold"))

dev.off()




##############################################################################
# # HISTOGRAMS PER BIOMARKER
# for (i in 7:2584){
#   hist(log2(as.numeric((data0[,i]))+1), main=colnames(data0)[i],xlab="number of reads in CPM (log2)",col="lightblue")
# }
# dev.off()
# 
# pdf("distribution of reads by patient all.pdf")
# 
#
# # PIE - READS OF BIOMARKERS WITHIN EACH PATIENT
# for (i in 1:188){
#   pie(as.numeric(data0[i,7:2584]),labels=colnames(data0)[7:2584],main=paste(as.character(data0[i,1]),"from", as.character(data0$group[i])))
# }
# dev.off()
# 
# # FILTER BIOMARKERS WITH MEDIAN IN AT LEAST ONE GROUP > 5 CPM
# # COMPUTE MEDIAN
# percentile=NULL
# for( i in 7:2584){
#   fifty=by(data0[,i],data0$group,median)
#   percentile=rbind(percentile,c(colnames(data0)[i],fifty))
# }
# e=percentile[which(as.numeric(percentile[,2])>5 | as.numeric(percentile[,3])>5 | as.numeric(percentile[,4])>5 | as.numeric(percentile[,5])>5
#                    | as.numeric(percentile[,6])>5 | as.numeric(percentile[,7])>5 
# ),1]
# 
# newdata=(data0[,c(1:6)])
# a=colnames(data0)
# 
# for( j in 1:length(e)) {
#   newdata=cbind(newdata,data0[,which(a==e[j])])
# }
# colnames(newdata)[-c(1:6)]=e
# data1=as.data.frame(newdata)
# 
# Group1=data1[which(data1$group=="Group 1"),]
# Group2=data1[which(data1$group=="Group 2"),]
# Group3a=data1[which(data1$group=="Group 3a"),]
# Group3b=data1[which(data1$group=="Group 3b"),]
# Group4a=data1[which(data1$group=="Group 4a"),]
# Group4b=data1[which(data1$group=="Group 4b"),]
# 
# ######Check the normality and the equality of the variance
# 
# n=ncol(Group1)
# normality(7,1580,Group1,"Normality of the data on Group1")
# normality(7,1580,Group2,"Normality of the data on Group2")
# normality(7,1580,Group3a,"Normality of the data on Group3a")
# normality(7,1580,Group3b,"Normality of the data on Group3b")
# normality(7,1580,Group4a,"Normality of the data on Group4a")
# normality(7,1580,Group4b,"Normality of the data on Group4b")
# 
# library(lawstat)
# test=NULL
# for(i in 7:1580){
#   test=rbind(test,c(colnames(data1)[i],levene.test(data1[,i],data1$group)$p.value))
# }
# write.csv(test,"equality of the variance.csv",row.names=F)
# 
# ################ANOVA
# 
# ANOVA_Non_param(7,1580,data1,data1$group,"Non parametric ANOVA for the group")
