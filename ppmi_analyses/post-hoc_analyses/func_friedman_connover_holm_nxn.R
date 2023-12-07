# Title: func_friedman_connover_holm_nxn.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains a function to perform Friedman test with Holm Post-hoc analysis for multiple pairwise comparisons (n x n comparisons)


friedman_connover_holm_nxn <- function(groups=NULL, blocks=NULL, score=NULL)  {
  testdata <- data.frame(
    "groups" = as.factor(groups),
    "blocks" = as.factor(blocks),
    "score" = as.double(score)
  )
  
  print("Friedman test:")
  friedman <- friedman.test(y = testdata$score,
                            groups = testdata$groups,
                            blocks = testdata$blocks)
  print(friedman)
  
  if (friedman$p.value <= 0.05) {
    print(paste("Friedman test results in significant differences among methods with respect to scores. P-value is", friedman$p.value))
    
    pvals <- frdAllPairsConoverTest(y = testdata$score,
                                    groups = testdata$groups,
                                    blocks = testdata$blocks,
                                    p.adjust.method = "holm")$p.value
    
    COMP_PERFORMANCE = list()
    for (i in colnames(pvals)){
      COMP_PERFORMANCE[[i]] = names(which(pvals[, i] < 0.05))
    }
    for (i in rownames(pvals)) {
      COMP_PERFORMANCE[[i]] = unique(c(COMP_PERFORMANCE[[i]], names(which(pvals[i, ] < 0.05))))
    }
    
    # generate ranks -----------
    avgranks = rep(0,length(levels(testdata$groups)))
    names(avgranks) <- levels(testdata$groups)
    for (i in levels(testdata$blocks)) {
      data2 <- 1 - testdata[testdata$blocks == i,]$score
      names(data2) <- testdata[testdata$blocks == i,]$groups
      ranks = rank(data2)
      for (j in levels(testdata$groups)) {
        avgranks[j] = avgranks[j] + ranks[j]
      }
    } 
    avgranks <- avgranks/length(levels(testdata$blocks))
    print(avgranks)
    
    # best and worst cases -----------
    best <- names(which.min(avgranks))
    if (length(COMP_PERFORMANCE[[best]])) {
      print(paste("Conover's all-pairs comparisons tests of Friedman-type ranked data with holm p-adjusted values,", best, 
                  "performs significantly (p = 0.05) better than", paste(COMP_PERFORMANCE[[best]], collapse=", ")))
    } else {
      print(paste("Conover's all-pairs comparisons tests of Friedman-type ranked data with holm p-adjusted values,", best, 
                  "didn't perform significantly (p = 0.05) better than any other of its category"))
    }
    worst <- names(which.max(avgranks))
    if (length(COMP_PERFORMANCE[[worst]])) {
      print(paste("Conover's all-pairs comparisons tests of Friedman-type ranked data with holm p-adjusted values,", worst, 
                  "performs significantly (p = 0.05) worse than", paste(COMP_PERFORMANCE[[worst]], collapse=", ")))
    } else {
      print(paste("Conover's all-pairs comparisons tests of Friedman-type ranked data with holm p-adjusted values,", worst, 
                  "didn't perform significantly (p = 0.05) worse than any other of its category"))
    }
    
  }
  else {
    print(paste("Friedman test does not show strong evidence agaist the null that groups are equivalent with respect to scores. P-value is", friedman$p.value))
  }

}
