library(rentrez)
getwd()
setwd("/home/jarek/extracting_sequences")

genomes_prok <- read.csv("genomes_proks.csv")
genomes_prok
head(genomes_prok$GenBank.FTP)
ftp_levels <- levels(genomes_prok$GenBank.FTP)
maxlev <- nlevels(genomes_prok$GenBank.FTP)
maxlev

wektor <- seq(1,4719)

setwd("/home/jarek/extracting_sequences/Genomes")

for (i in wektor){
  x <- ftp_levels[[i]];
  z <- strsplit(x, "/");
  OMG <- paste0(x,"/", z[[1]][10],"_genomic.gbff.gz");
  cat(i);
  # tt <- list.files()
  # if(is.element(OMG, tt) == T){cat("skipping file")}else{download.file(url = OMG, destfile = z[[1]][10])}
  cat(download.file(url = OMG, destfile = paste0(z[[1]][10], "_genomic.gbff.gz")))
  wektor <- setdiff(wektor, i)
}
