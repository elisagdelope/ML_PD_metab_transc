# Title: func_friedman.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains a function to perform Friedman test (friedman); 
# a function to perform post-hoc analysis with Holm correction for many-to-one comparisons (friedman_holm_nx1);
# and a function to perform post-hoc analysis with Bergmann-Hommel correction for multiple pairwise comparisons (friedman_BergmannHommel_nxn)

library(PMCMRplus)
library(ggplot2)
library(scmamp)

friedman <- function(groups=NULL, blocks=NULL, score=NULL)  {
  testdata <- data.frame(
    "groups" = as.factor(groups),
    "blocks" = as.factor(blocks),
    "score" = as.double(score)
  )
  
  frd <- friedman.test(y = testdata$score,
                       groups = testdata$groups,
                       blocks = testdata$blocks)
  print(frd)
  
  if (frd$p.value <= 0.05) {
    print(paste("Friedman test results in significant differences among methods with respect to scores. P-value is", frd$p.value))
  }else {
    print(paste("Friedman test does not show strong evidence agaist the null that groups are equivalent with respect to scores. P-value is", frd$p.value))
  }
  return(frd)
}



friedman_holm_nx1 <- function(groups=NULL, blocks=NULL, score=NULL, col_control=1)  {
  testdata <- data.frame(
    "groups" = as.factor(groups),
    "blocks" = as.factor(blocks),
    "score" = as.double(score)
  )
  testdata <- testdata %>% 
    pivot_wider(names_from=groups, values_from = score)
  
  rm <- colMeans(scmamp::rankMatrix(testdata[,-1]))
  print(rm)
  # many-to-one (nx1)
  pv <- friedmanPost(as.matrix(testdata[,-1]), control=col_control)
  pv.adj <- p.adjust(pv, method = "holm", )
  
  df <- as.data.frame(rbind(pv, pv.adj))
  rownames(df) <- c("PV", "PAV-holm")
  
  print(df)
  return(list(test=df, ranks=rm))
}



friedman_BergmannHommel_nxn <- function(groups=NULL, blocks=NULL, score=NULL)  {
  testdata <- data.frame(
    "groups" = as.factor(groups),
    "blocks" = as.factor(blocks),
    "score" = as.double(score)
  )
  testdata <- testdata %>% 
    pivot_wider(names_from=groups, values_from = score)
  # all pairwise comparisons (nxn)
  pv.matrix <- friedmanPost(as.matrix(testdata[,-1]), control=NULL)
  pv.adj.matrix <- adjustBergmannHommel(pv.matrix)
  r.means <- colMeans(scmamp::rankMatrix(as.matrix(testdata[,-1])))
  print(r.means)
  print(pv.adj.matrix)
  
  #plot
  plt <- plotPvalues(pvalue.matrix=pv.adj.matrix, alg.order=order(r.means, decreasing=TRUE), font.size = 5)
  plt <- plt + 
    labs(title="Friedman test with adjusted PV (Bergmann-Hommel)") + 
    scale_fill_gradientn("Adj. PV" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
    labs(x = "Feature type", y = "Feature type") + theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),    # Set x-axis label size
      axis.title.y = element_text(size = 14),     # Set y-axis label size
      plot.title = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) 
  return(list(test=pv.adj.matrix, plot = plt, ranks=r.means))
}
