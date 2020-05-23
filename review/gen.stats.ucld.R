library("HDInterval")

tmrca.all <- read.csv("seqs01/all-true-tmrca.tab", header = FALSE)[[1]]
ucld.dir <- "ucld/"
last.sample <- 2020

logs.dir <- paste(ucld.dir,"logs/",sep='')


# ucld.data <- read.table("<logs.dir>/101.ucld.1.log", sep="\t", header=TRUE)
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



df.rates <- data.frame(matrix(ncol=6,nrow=0))
df.tmrca <- data.frame(matrix(ncol=7,nrow=0))

burnin <- 0.1

for(tree in seq(101,200)) {
  # row <- c()
  tmrca.true <- tmrca.all[tree-100]
  ucld.file <-  paste(logs.dir,tree,".ucld.1.log",sep='')
  # ucld samples
  log.data <- read.table(ucld.file,sep="\t", header=TRUE)
  n <- nrow(log.data)
  s <- as.integer(burnin*n)
  log.data <- log.data[s:n, ]
  ucld.rates <- log.data["ucldMean"][[1]]
  urate.mean <- mean(ucld.rates)
  urate.sd <- sd(ucld.rates)
  hd95 <- hdi(ucld.rates, credMass=0.95)
  q95 <- quantile(ucld.rates, c(0.05,0.95))
  rates.row <- c(urate.mean, urate.sd, hd95["lower"][[1]],hd95["upper"][[1]],q95[[1]],q95[[2]])
  df.rates <- rbind(df.rates, rates.row)  
  ucld.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  utmrca.mean <- mean(ucld.tmrcas)
  utmrca.sd <- sd(ucld.tmrcas)
  hd95 <- hdi(ucld.tmrcas, credMass=0.95)
  q95 <- quantile(ucld.tmrcas, c(0.05,0.95))  
  tmrca.row <- c(tmrca.true, utmrca.mean, utmrca.sd, hd95["lower"][[1]],hd95["upper"][[1]],q95[[1]],q95[[2]] )
  df.tmrca <- rbind(df.tmrca, tmrca.row)
  print(rates.row)
}

colnames(df.rates) <- c("urate.mean","urate.sd","urate.l","urate.u","urate.q5","urate.q95")
colnames(df.tmrca) <- c("tmrca.true", "utmrca.mean", "utmrca.sd", "utmrca.hd95.l","utmrca.hd95.u","utmrca.q5","utmrca.q95" )


write.csv(df.rates, file=paste(ucld.dir,"stats_rates.csv",sep=''), row.names=FALSE)

write.csv(df.tmrca, file=paste(ucld.dir,"stats_tmrca.csv",sep=''), row.names=FALSE)


