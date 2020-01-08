# Assumption: all stat files have 100 entries corresponding to trees 1..100

ucld.rates <- read.csv("../test5.ucld/stats_rates.csv", header = TRUE)
arc6.rates <-  read.csv("../test6/stats_rates.csv", header = TRUE)
arc7.rates <-  read.csv("../test7/stats_rates.csv", header = TRUE)

df.rates <- data.frame(matrix(ncol=0,nrow=100))
df.rates$tree <- seq(1,100)
df.rates$ucld.rate <- ucld.rates$urate.mean
df.rates$arc6.rate <- arc6.rates$arate.mean
df.rates$arc7.rate <- arc7.rates$arate.mean
df.rates$arc6.omega <- arc6.rates$aomega.mean
df.rates$arc7.omega <- arc7.rates$aomega.mean


write.csv(df.rates, file="mean_rates.csv", row.names=FALSE)

tmrca.all <- read.csv("../tmrca.tab", header = FALSE)[[1]]
ucld.tmrca <- read.csv("../test5.ucld/stats_tmrca.csv", header = TRUE)
arc6.tmrca <-  read.csv("../test6/stats_tmrca.csv", header = TRUE)
arc7.tmrca <-  read.csv("../test7/stats_tmrca.csv", header = TRUE)

df.tmrca <- data.frame(matrix(ncol=0,nrow=100))
df.tmrca$tree <- seq(1,100)
df.tmrca$truev <- tmrca.all
df.tmrca$ucld <- ucld.tmrca$utmrca.mean
df.tmrca$arc6 <- arc6.tmrca$atmrca.mean
df.tmrca$arc7 <- arc7.tmrca$atmrca.mean

write.csv(df.tmrca, file="mean_tmrca.csv", row.names=FALSE)