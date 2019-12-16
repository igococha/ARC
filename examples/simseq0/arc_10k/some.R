library(ape)

x <- read.nexus.data("1.nex")

names(x)

paste(x['99_2009.8'][[1]], collapse='',sep='' )

paste(x[names(x)[[2]]][[1]] , collapse='',sep='' )


