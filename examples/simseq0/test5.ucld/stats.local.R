
tmrca.all <- read.csv("../tmrca.tab", header = FALSE)[[1]]
ucld.dir <- "ucld_logs/"
last.sample <- 2010

# ucld.data <- read.table("ucld_hpc_logs/10.ucld.5.log", sep="\t", header=TRUE)
# > names(ucld.data)
#  [1] "Sample"                      "posterior"                  
#  [3] "likelihood"                  "prior"                      
#  [5] "treeLikelihood"              "TreeHeight"                 
#  [7] "ucldMean"                    "ucldStdev"                  
#  [9] "kappa"                       "mutationRate"               
# [11] "gammaShape"                  "popSize"                    
# [13] "CoalescentConstant"          "rate.mean"                  
# [15] "rate.variance"               "rate.coefficientOfVariation"
# [17] "freqParameter.1"             "freqParameter.2"            
# [19] "freqParameter.3"             "freqParameter.4"     

# arc.data <- read.table("arc_hpc_logs/10.arc.5.log", sep="\t", header=TRUE)
# > names(arc.data)
#  [1] "Sample"             "posterior"          "likelihood"        
#  [4] "prior"              "treeLikelihood"     "TreeHeight"        
#  [7] "rateOmega"          "rateMean"           "kappa"             
# [10] "mutationRate"       "gammaShape"         "popSize"           
# [13] "CoalescentConstant" "freqParameter.1"    "freqParameter.2"   
# [16] "freqParameter.3"    "freqParameter.4"  


df.rates <- data.frame(matrix(ncol=6,nrow=0))
df.tmrca <- data.frame(matrix(ncol=5,nrow=0))

for(tree in seq(1,100)) {
  # row <- c()
  tmrca.true <- tmrca.all[tree]
  ucld.file <-  paste(ucld.dir,tree,".ucld.5.log",sep='')
  # ucld samples
  log.data <- read.table(ucld.file,sep="\t", header=TRUE)
  n <- nrow(log.data)
  s <- as.integer(0.1*n)
  log.data <- log.data[s:n, ]
  ucld.rates <- log.data["ucldMean"][[1]]
  urate.mean <- mean(ucld.rates)
  urate.sd <- sd(ucld.rates)
  ucld.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  utmrca.mean <- mean(ucld.tmrcas)
  utmrca.sd <- sd(ucld.tmrcas)
  rates.row <- c(urate.mean, urate.sd)
  df.rates <- rbind(df.rates, rates.row)
  tmrca.row <- c(tmrca.true, utmrca.mean, utmrca.sd)
  df.tmrca <- rbind(df.tmrca, tmrca.row)
  print(rates.row)
}

colnames(df.rates) <- c("urate.mean","urate.sd")
colnames(df.tmrca) <- c("tmrca.true", "utmrca.mean", "utmrca.sd")


write.csv(df.rates, file="stats_rates.csv", row.names=FALSE)

write.csv(df.tmrca, file="stats_tmrca.csv", row.names=FALSE)