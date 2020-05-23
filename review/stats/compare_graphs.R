# Assumption: all stat files have 100 entries corresponding to trees 101..200

library(ggplot2)
library(reshape)

ucld.rates <- read.csv("../ucld/stats_rates.csv", header = TRUE)
arc.rates <-  read.csv("../arc/stats_rates.csv", header = TRUE)


df.rates <- data.frame(matrix(ncol=0,nrow=100))
df.rates$tree <- seq(101,200)
df.rates$truev <- 0.0005
df.rates$ucld <- ucld.rates$urate.mean
df.rates$arc <- arc.rates$arate.mean


df.rates <- melt(df.rates, id=c("tree"))
names(df.rates) <- c("tree","method","mean.rate")
df.rates$tree <- as.numeric(as.character(df.rates$tree))


png("rates.alltrees.png")
ggplot(df.rates, aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU - all trees")
dev.off()

png("rates.trees_101_125.png")
ggplot(df.rates[df.rates$tree < 126,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU : 101-125")
dev.off()

png("rates.trees_126_150.png")
ggplot(df.rates[df.rates$tree > 125 & df.rates$tree < 151,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU : 126-150")
dev.off()

png("rates.trees_151_175.png")
ggplot(df.rates[df.rates$tree > 150 & df.rates$tree < 176,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU: 151-175")
dev.off()

png("rates.trees_176_200.png")
ggplot(df.rates[df.rates$tree > 175 ,], aes(x = mean.rate, y = tree, colour = method)) +
  geom_point() + ggtitle("mean MU: 176-200")
dev.off()

### tmrca

tmrca.all <- read.csv("../seqs01/all-true-tmrca.tab", header = FALSE)[[1]]
ucld.tmrca <- read.csv("../ucld/stats_tmrca.csv", header = TRUE)
arc.tmrca <-  read.csv("../arc/stats_tmrca.csv", header = TRUE)

df.tmrca <- data.frame(matrix(ncol=0,nrow=100))
df.tmrca$tree <- seq(101,200)
df.tmrca$truev <- tmrca.all
df.tmrca$ucld <- ucld.tmrca$utmrca.mean
df.tmrca$arc <- arc.tmrca$atmrca.mean

df.tmrca <- melt(df.tmrca, id=c("tree"))
names(df.tmrca) <- c("tree","method","tmrca")

png("tmrca.alltrees.png")
ggplot(df.tmrca, aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA - all trees")
dev.off()

png("tmrca.trees_101_125.png")
ggplot(df.tmrca[df.tmrca$tree < 126,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 101-125")
dev.off()

png("tmrca.trees_126_150.png")
ggplot(df.tmrca[df.tmrca$tree > 125 & df.tmrca$tree < 151,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 126-150")
dev.off()

png("tmrca.trees_151_175.png")
ggplot(df.tmrca[df.tmrca$tree > 150 & df.tmrca$tree < 176,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 151-175")
dev.off()

png("tmrca.trees_176_200.png")
ggplot(df.tmrca[df.tmrca$tree > 175 ,], aes(x = tmrca, y = tree, colour = method)) +
  geom_point() + ggtitle("TMRCA: 176-200")
dev.off()
