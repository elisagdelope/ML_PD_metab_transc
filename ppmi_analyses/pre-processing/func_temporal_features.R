# Title: func_temporal_features.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains functions that transform a t number of features measured over several timepoints (t dimensions) into temporal features of 1 dimension (ratio, sd, difference,...)

library(dplyr)
library(tidyr)


# standard deviation over first visits (needs >=3 visits)
exprs_sd <- function(expression, varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "SAMPLE", values_to = "EXPRS") %>%
    separate(SAMPLE, c("PATIENT_ID", "VISIT"), sep = '\\.') %>%
    group_by(!!as.symbol(varname), PATIENT_ID) %>% 
    summarise(EXPRS = sd(EXPRS)) %>% # sd
    pivot_wider(names_from = all_of(varname), values_from = EXPRS) 
  return(expression)
}

# lm (exprs ~ time) (needs >=3 visits)
exprs_lm <- function(expression, varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "SAMPLE", values_to = "EXPRS") %>%
    separate(SAMPLE, c("PATIENT_ID", "VISIT"), sep = '\\.') %>%
    group_by(!!as.symbol(varname), PATIENT_ID) %>% 
    mutate(TIMEPOINT = as.integer(case_when(
      as.character(VISIT) == "BL" ~ 1,
      as.character(VISIT) == "V04" ~ 2,
      as.character(VISIT) == "V06" ~ 3,
      as.character(VISIT) == "V08" ~ 4))) %>% 
    summarise(EXPRS = lm(EXPRS ~ TIMEPOINT)$coefficients[[2]]) %>% # slope coefficient of lm (exprs ~ time)
    pivot_wider(names_from = all_of(varname), values_from = EXPRS) 
  return(expression)
}

# lm (exprs ~ lag) (needs >=4 visits)
exprs_lm_lag <- function(expression, varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "SAMPLE", values_to = "EXPRS") %>%
    separate(SAMPLE, c("PATIENT_ID", "VISIT"), sep = '\\.') %>%
    group_by(!!as.symbol(varname), PATIENT_ID) %>% 
    mutate(LAG1 = lag(EXPRS)) %>% 
    drop_na() %>% # remove first timepoint (lag1 = NA)
    summarise(EXPRS = lm(EXPRS ~ LAG1)$coefficients[[2]]) %>% # slope coefficient of lm (exprs ~ lag1)
    pivot_wider(names_from = all_of(features_varname), values_from = EXPRS)
  return(expression)
}




# computes ratio of expression on t2 over t1: t2/t1 from expression matrix
exprs_ratio_t2t1 <- function(expression, visits_list, features, features_varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "exprs") %>%
    separate(sample, c("patient", "visit"), sep = '\\.') %>%
    filter((visit %in% visits_list) & (get(features_varname) %in% features)) %>%
    group_by(!!as.symbol(features_varname), patient) %>% 
    filter(n() == length(visits_list)) %>%   # remove patients that only have data for 1 visit
    summarise(exprs = ifelse(exprs[visit == visits_list[1]] != 0, exprs[visit == visits_list[2]] / exprs[visit == visits_list[1]], exprs[visit == visits_list[2]] / 0.00001), .groups = "keep") %>% # ratio
    pivot_wider(names_from = all_of(features_varname), values_from = exprs)
  return(expression)
}


# standard deviation over first visits (needs 3 visits)
exprs_sd_t3 <- function(expression, visits_list, features, features_varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "exprs") %>%
    separate(sample, c("patient", "visit"), sep = '\\.') %>%
    filter((visit %in% visits_list) & (get(features_varname) %in% features)) %>%
    group_by(!!as.symbol(features_varname), patient) %>% 
    filter(n() == length(visits_list)) %>%   # remove patients that only have data for 1 visit
    summarise(exprs = sd(exprs)) %>% # sd
    pivot_wider(names_from = all_of(features_varname), values_from = exprs)
  return(expression)
}



# absolute difference in relation to median absolute difference across all genes -> GENE OR PATIENT??
exprs_sd_t3 <- function(expression, visits_list, features, features_varname) {
  expression <- expression %>%
    pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "exprs") %>%
    separate(sample, c("patient", "visit"), sep = '\\.') %>%
    filter((visit %in% visits_list) & (get(features_varname) %in% features)) %>%
    group_by(!!as.symbol(features_varname), patient) %>% 
    filter(n() == length(visits_list)) %>%   # remove patients that only have data for 1 visit
    summarise(exprs = abs(exprs[visit == visits_list[2]] - exprs[visit == visits_list[1]])) %>% # absolute diff
    pivot_wider(names_from = features_varname, values_from = exprs)
  return(expression)
}