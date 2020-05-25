library("HDInterval")

tmrca.all <- read.csv("seqs01/all-true-tmrca.tab", header = FALSE)[[1]]
arc.dir <- "arc/"
# ucld.dir <- "ucld/logs/"
last.sample <- 2020

arc.logs.dir <- paste(arc.dir,"logs/",sep='')


# arc.data <- read.table("<logs.dir>/101.arc.1.log", sep="\t", header=TRUE)
# > names(arc.data)
#  [1] "Sample"             "posterior"          "likelihood"        
#  [4] "prior"              "treeLikelihood"     "TreeHeight"        
#  [7] "rateOmega"          "rateMean"           "kappa"             
# [10] "mutationRate"       "gammaShape"         "popSize"           
# [13] "CoalescentConstant" "freqParameter.1"    "freqParameter.2"   
# [16] "freqParameter.3"    "freqParameter.4"  


df.rates <- data.frame(matrix(ncol=8,nrow=0))
df.tmrca <- data.frame(matrix(ncol=7,nrow=0))

burnin <- 0.1

for(tree in seq(101,200)) {
  # row <- c()
  tmrca.true <- tmrca.all[tree-100]
  # ucld samples
  #ucld.file <-  paste(ucld.dir,tree,".ucld.1.log",sep='')  
  #log.data <- read.table(ucld.file,sep="\t", header=TRUE)
  #n <- nrow(log.data)
  #s <- as.integer(burnin*n)
  #log.data <- log.data[s:n, ]
  #ucld.rates <- log.data["ucldMean"][[1]]
  #urate.mean <- mean(ucld.rates)
  #urate.sd <- sd(ucld.rates)
  #ucld.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  #utmrca.mean <- mean(ucld.tmrcas)
  #utmrca.sd <- sd(ucld.tmrcas)

  # arc samples
  arc.file <- paste(arc.logs.dir,tree,".arc.1.log",sep='')
  log.data <- read.table(arc.file,sep="\t", header=TRUE)
  n <- nrow(log.data)
  s <- as.integer(burnin*n)
  log.data <- log.data[s:n, ]
  arc.rates <- log.data["rateMean"][[1]]
  arate.mean <- mean(arc.rates)
  arate.sd <-  sd(arc.rates)
  hd95 <- hdi(arc.rates, credMass=0.95)
  q95 <- quantile(arc.rates, c(0.05,0.95))
  arc.omegas <- log.data["rateOmega"][[1]]
  aomega.mean <- mean(arc.omegas)
  aomega.sd <- sd(arc.omegas)
  aomega.hd95 <- hdi(arc.omegas, credMass=0.95)
  aomega.q95 <- quantile(arc.omegas, c(0.05,0.95))
  # rates.row <- c(urate.mean, urate.sd, arate.mean, arate.sd, aomega.mean, aomega.sd)
  rates.row <- c(arate.mean, arate.sd, hd95["lower"][[1]],hd95["upper"][[1]],q95[[1]],q95[[2]], aomega.mean, aomega.sd,aomega.hd95["lower"][[1]],aomega.hd95["upper"][[1]],aomega.q95[[1]],aomega.q95[[2]])
  df.rates <- rbind(df.rates, rates.row)

  arc.tmrcas <- last.sample - log.data["TreeHeight"][[1]]
  atmrca.mean <- mean(arc.tmrcas)
  atmrca.sd <- sd(arc.tmrcas)
  hd95 <- hdi(arc.tmrcas, credMass=0.95)
  q95 <- quantile(arc.tmrcas, c(0.05,0.95))  

  tmrca.row <- c(tmrca.true, atmrca.mean, atmrca.sd, hd95["lower"][[1]],hd95["upper"][[1]],q95[[1]],q95[[2]])
  df.tmrca <- rbind(df.tmrca, tmrca.row)
  print(rates.row)
}

colnames(df.rates) <- c("arate.mean","arate.sd","arate.l","arate.u","arate.q5","arate.q95" ,"aomega.mean","aomega.sd","aomega.l","aomega.u","aomega.q5","aomega.q95")
colnames(df.tmrca) <- c("tmrca.true","atmrca.mean","atmrca.sd","atmrca.hd95.l","atmrca.hd95.u","atmrca.q5","atmrca.q95" )

write.csv(df.rates, file=paste(arc.dir,"stats_rates.csv",sep=''), row.names=FALSE)
write.csv(df.tmrca, file=paste(arc.dir,"stats_tmrca.csv",sep=''), row.names=FALSE)

