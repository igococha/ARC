
tmrca.all <- read.csv("tmrca.tab", header = FALSE)[[1]]
arc.dir <- "arc_hpc_logs/"
ucld.dir <- "ucld_hpc_logs/"
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


df.rates <- data.frame(matrix(ncol=4,nrow=0))
df.tmrca <- data.frame(matrix(ncol=5,nrow=0))


for(tree in seq(1,100)) {
  # row <- c()
  tmrca.true <- tmrca.all[tree]
  ucld.file <-  paste(ucld.dir,tree,".ucld.5.log",sep='')
  # ucld samples
  log.data <- read.table(ucld.file,sep="\t", header=TRUE)
  ucld.rates <- log.data["ucldMean"][[1]]
  urate.mean <- mean(ucld.rates)
  urate.sd <- sd(ucld.rates)
  ucld.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  utmrca.mean <- mean(ucld.tmrcas)
  utmrca.sd <- sd(ucld.tmrcas)
  # arc samples
  arc.file <- paste(arc.dir,tree,".arc.5.log",sep='')
  log.data <- read.table(arc.file,sep="\t", header=TRUE)
  arc.rates <- log.data["rateMean"][[1]]
  arate.mean <- mean(arc.rates)
  arate.sd <-  sd(arc.rates)
  arc.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  atmrca.mean <- mean(arc.tmrcas)
  atmrca.sd <- sd(arc.tmrcas)
  rates.row <- c(urate.mean, urate.sd, arate.mean, arate.sd)
  df.rates <- rbind(df.rates, rates.row)
  tmrca.row <- c(tmrca.true, utmrca.mean, utmrca.sd,  atmrca.mean, atmrca.sd)
  df.tmrca <- rbind(df.tmrca, tmrca.row)
  print(rates.row)
}

colnames(df.rates) <- c("urate.mean","urate.sd","arate.mean","arate.sd")
colnames(df.tmrca) <- c("tmrca.true", "utmrca.mean", "utmrca.sd",
                        "atmrca.mean","atmrca.sd")


write.csv(df.rates, file="stats_rates.csv", row.names=FALSE)

write.csv(df.tmrca, file="stats_tmrca.csv", row.names=FALSE)