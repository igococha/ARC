# Assumption: all stat files have 100 entries corresponding to trees 101..200

ucld.rates <- read.csv("../ucld/stats_rates.csv", header = TRUE)
arc.rates <-  read.csv("../arc/stats_rates.csv", header = TRUE)


df.rates <- data.frame(matrix(ncol=0,nrow=100))
df.rates$tree <- seq(101,200)
df.rates$ucld.rate <- ucld.rates$urate.mean
df.rates$arc.rate <- arc.rates$arate.mean
df.rates$arc.omega <- arc.rates$aomega.mean

write.csv(df.rates, file="mean_rates.csv", row.names=FALSE)


tmrca.all <- read.csv("../seqs01/all-true-tmrca.tab", header = FALSE)[[1]]
ucld.tmrca <- read.csv("../ucld/stats_tmrca.csv", header = TRUE)
arc.tmrca <-  read.csv("../arc/stats_tmrca.csv", header = TRUE)

df.tmrca <- data.frame(matrix(ncol=0,nrow=100))
df.tmrca$tree <- seq(101,200)
df.tmrca$truev <- tmrca.all
df.tmrca$ucld <- ucld.tmrca$utmrca.mean
df.tmrca$arc <- arc.tmrca$atmrca.mean

write.csv(df.tmrca, file="mean_tmrca.csv", row.names=FALSE)