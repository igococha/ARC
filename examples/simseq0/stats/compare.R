# Assumption: all stat files have 100 entries corresponding to trees 1..100

library(ggplot2)
library(reshape)

ucld.rates <- read.csv("../test5.ucld/stats_rates.csv", header = TRUE)
arc6.rates <-  read.csv("../test6/stats_rates.csv", header = TRUE)
arc7.rates <-  read.csv("../test7/stats_rates.csv", header = TRUE)

df.rates <- data.frame(matrix(ncol=0,nrow=100))
df.rates$tree <- seq(1,100)
df.rates$truev <- 0.005
df.rates$ucld <- ucld.rates$urate.mean
df.rates$arc6 <- arc6.rates$arate.mean
df.rates$arc7 <- arc7.rates$arate.mean

df.rates <- melt(df.rates, id=c("tree"))
names(df.rates) <- c("tree","method","mean.rate")
df.rates$tree <- as.numeric(as.character(df.rates$tree))

png("rates.alltrees.png")
ggplot(df.rates, aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU - all trees")
dev.off()

png("rates.trees_0_25.png")
ggplot(df.rates[df.rates$tree < 26,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU : 1-25")
dev.off()

png("rates.trees_26_50.png")
ggplot(df.rates[df.rates$tree > 25 & df.rates$tree < 51,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU : 26-50")
dev.off()

png("rates.trees_51_75.png")
ggplot(df.rates[df.rates$tree > 50 & df.rates$tree < 76,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU: 51-75")
dev.off()

png("rates.trees_76_100.png")
ggplot(df.rates[df.rates$tree > 75 ,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU: 76-100")
dev.off()

### tmrca

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

df.tmrca <- melt(df.tmrca, id=c("tree"))
names(df.tmrca) <- c("tree","method","tmrca")

png("tmrca.alltrees.png")
ggplot(df.tmrca, aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA - all trees")
dev.off()

png("tmrca.trees_0_25.png")
ggplot(df.tmrca[df.tmrca$tree < 26,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 1-25")
dev.off()

png("tmrca.trees_26_50.png")
ggplot(df.tmrca[df.tmrca$tree > 25 & df.tmrca$tree < 51,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 26-50")
dev.off()

png("tmrca.trees_51_75.png")
ggplot(df.tmrca[df.tmrca$tree > 50 & df.tmrca$tree < 76,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 51-75")
dev.off()

png("tmrca.trees_76_100.png")
ggplot(df.tmrca[df.tmrca$tree > 75 ,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 76-100")
dev.off()
