library(ape)

seq.file <- "1.nex"

seq.nexus <- read.nexus.data(seq.file)

tip.labels <- names(seq.nexus)

tip.dates <- c()
for(tip.label in tip.labels) {
  tip.date <- strsplit(tip.label,'_')[[1]][[2]]
  tip.dates <- c(tip.dates, tip.date)
}
tip.dates <- sort(tip.dates)

first.sample <- tip.dates[1]
last.sample <- tip.dates[length(tip.dates)]

print(paste("First sample:",first.sample))
print(paste("Last sample:",last.sample))