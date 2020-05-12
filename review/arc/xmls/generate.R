# Script for generating beast xml files from nexus sequences using an xml template file.

library(ape)


### Parameters
# sequences directory
seq.dir <- "../../seqs01/"
template.xml <- "arc.template"
type <- "arc"

xml.version <- 1

# --numsites--
# --sequences--
# --logfile--

genXMLSequence <- function(seq.nexus) {
  tip.labels <- names(seq.nexus)
  seq.xml <- ""
  date.trait <- c()
  for(tip.label in tip.labels) {
    tip.date <- strsplit(tip.label,'_')[[1]][[2]]
    seq <- paste(seq.nexus[tip.label][[1]], collapse='', sep='')
    seq.len <- nchar(seq)
    seq.xml <- paste(seq.xml,"<sequence id='", tip.label,"' spec='Sequence'",sep='')
    seq.xml <- paste(seq.xml, " taxon='",tip.label,"' totalcount='4' value=",sep='')
    seq.xml <- paste(seq.xml,"'",seq,"'/>\n",sep='')
    date.trait <- c(date.trait,paste(tip.label,"=",tip.date,sep=''))
  }
  date.trait <- paste(date.trait,collapse=",")
  return(c(seq.len,seq.xml,date.trait))
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

read.my.dna <- function(ff,start=2000) {
  con = file(ff, "r")
  l <- trim(readLines(con,n=1))
  lsplit <- strsplit(l,' ')
  ntax <- as.numeric(lsplit[[1]][[1]])
  nchar <- as.numeric(lsplit[[1]][[2]])
  i <- 0
  res <- c()
  while ( TRUE ) {
    l = readLines(con, n = 1)
    if ( length(l) == 0 ) {
      break
    }
    l <- trim(l)
    lsplit <- strsplit(l,"\\s+")
    label <- lsplit[[1]][[1]]
    seq <- lsplit[[1]][[2]]
    seq.date <- start+(as.numeric(label)-1)/10
    i <- i+1
    label <- paste(label,'_',seq.date,sep='')
    res[label] <- seq
  }
  close(con)
  return(res)
}

seq.files <- list.files(path=seq.dir,pattern="*.nex", full.names=FALSE)

xmls <- c()
for(f in seq.files) {
  # seq.nexus <- read.nexus.data(paste(seq.dir,f,sep=''))
  seq.nexus <- read.my.dna(paste(seq.dir,f,sep=''))
  seq.info <- genXMLSequence(seq.nexus)
  seq.len <- seq.info[[1]]
  seq.xml <- seq.info[[2]]
  date.trait <- seq.info[[3]]
  tree.number <- strsplit(f,'\\.')[[1]][[1]]
  file.xml <- paste(tree.number,'.',type,'.',xml.version,'.xml',sep='')
  file.log <- paste(tree.number,'.',type,'.',xml.version,'.log',sep='')
  file.str <- readChar(template.xml,file.info(template.xml)$size)
  file.str <- gsub("--sequences--",seq.xml,file.str)
  file.str <- gsub("--numsites--",seq.len,file.str)
  file.str <- gsub("--logfile--",file.log,file.str)
  file.str <- gsub("--datetrait--",date.trait,file.str)
  sink(file.xml)
  cat(file.str)
  sink()
  xmls <- c(xmls,file.xml)
}

if (FALSE) {
  xmls <- c()
  for(f in seq.files) {
    tree.number <- strsplit(f,'\\.')[[1]][[1]]
    file.xml <- paste(tree.number,'.',type,'.',xml.version,'.xml',sep='')
    xmls <- c(xmls,file.xml)
  }
}

if (FALSE) {
  for(f in xmls) {
    system(paste("java -jar beastARC_2.jar",f))
  }
}

if (FALSE) {
  for(f in xmls) {
    system(paste("java -jar beastARC_2.jar -resume",f))
  }
}

# first check if dates are different
checkTipLabels <- function(seq.files) {
  all.tip.labels <- c()
  for(f in seq.files) {
    seq.nexus <- read.nexus.data(paste(seq.dir,f,sep=''))
    tip.labels <- names(seq.nexus)
    tip.labels <- sort(tip.labels)
    tip.labels <- paste(tip.labels,collapse='')
    all.tip.labels <- c(all.tip.labels, tip.labels)
  }
  # check first label set against rest
  bv <- (all.tip.labels != all.tip.labels[[1]])
  # which(bv)
}

# --sequences--
# --logfile--