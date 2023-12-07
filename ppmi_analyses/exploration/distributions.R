# fitting distributions
require("fitdistrplus")
set.seed(10)
n = 25
size = 27
prob = .4
data = rbinom(n, size = size, prob = prob)
fit = fitdist(data = data, dist="binom", 
              fix.arg=list(size = size), 
              start=list(prob = 0.1))

summary(fit)
plot(fit)
fit2 = fitdist(data, dist = "pois")
summary(fit2)
plot(fit2)

#-----
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PATHWAY = paste0("../", analysis_name , "/04-pathway_level") 
e_level = "GOCC"
st = "pathifier"
target = "DIAGNOSIS"
var_id = paste(e_level, "name", sep = "_")
EXPRESSION.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>% 
  column_to_rownames(var = var_id) 

AIC = list()
for (i in seq(from = 15, to=nrow(expression), by = 1)) {
  mydata = as.numeric(as.vector(expression[i,]))
  
  # binom
  f = fitdist(data = round(mydata*100)[1050:1200], dist="binom",
                      fix.arg=list(size = length(mydata[1050:1200])), 
                       start=list(prob = 0.1))
 # summary(f)
  plot(f)
  AIC[["binom"]] = f$aic
  print(paste("binom", AIC[["binom"]]))
  
  # nbinom
  f = fitdist(data = round(mydata*100)[1150:1200], dist="nbinom",
                       fix.arg=list(size = length(mydata[1150:1200])), 
                       start=list(prob = 0.1))
 # summary(f)
  plot(f)
  AIC[["nbinom"]] = f$aic
  print(paste("nbinom", AIC[["nbinom"]]))
  
  
  for (d in c("pois", "exp", "norm")) {
    fit_mydata2 = fitdist(data = round(mydata*100)[1150:1200], dist=d)
  #  summary(fit_mydata2)
    plot(fit_mydata2)
    AIC[[d]] = fit_mydata2$aic
    print(paste(d, AIC[[d]]))
  }
  
  
  mind = which(AIC == min(unlist(AIC)))
  print(paste("min AIC of distribution", mind) )
}







# test for binomial --------------------------------
for (db in c("GOBP", "GOCC", "CORUM")) {
  var_id = paste(db, "name", sep = "_")
  for (st in c("mean", "median", "sd", "pathifier")) { 
    file <- file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_"))
    expression <- vroom(file, col_types = cols(), delim = "\t") %>% # load data 
      column_to_rownames(var = var_id) 
    print(paste("analyzing", db, st))
    print(paste("# genes:", nrow(expression)))
    count_norm = 0 
    count_non_norm = 0
    for (i in seq(1:nrow(expression))) {
      mydata = as.numeric(as.vector(expression[i,]))
      fit_mydata = fitdist(data = round(mydata/10), dist="binom",
                           fix.arg=list(size = length(mydata)), 
                           start=list(prob = 0.1))
      
      if (s$p.value < 0.05) {
        count_non_norm = count_non_norm + 1
      } else {
        count_norm = count_norm + 1
      }
    }
    
    print(paste("# norm:", count_norm))
    print(paste("# non norm:", count_non_norm))
  }
}




# test for normality --------------------------------

for (db in c("GOBP", "GOCC", "CORUM")) {
  var_id = paste(db, "name", sep = "_")
  for (st in c("mean", "median", "sd", "pathifier")) { 
    file <- file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_"))
    expression <- vroom(file, col_types = cols(), delim = "\t") %>% # load data 
      column_to_rownames(var = var_id) 
    print(paste("analyzing", db, st))
    print(paste("# genes:", nrow(expression)))
    count_norm = 0 
    count_non_norm = 0
    for (i in seq(1:nrow(expression))) {
      mydata = as.numeric(as.vector(expression[i,]))
      s = shapiro.test(mydata)
#      print(paste(i, s$p.value))
      if (s$p.value < 0.05) {
        count_non_norm = count_non_norm + 1
      } else {
        count_norm = count_norm + 1
      }
    }
    
    print(paste("# norm:", count_norm))
    print(paste("# non norm:", count_non_norm))
  }
}

