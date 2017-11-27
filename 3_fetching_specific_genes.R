library(rentrez)
library(seqinr)

setwd("/home/jarek/extracting_sequences/Fragmented")
lista_plikow <- list.files("/home/jarek/extracting_sequences/Fragmented")
lista_plikow

for(i in length(lista_plikow)){
  plik_gbff <- lista_plikow[i];
  plik_fasta <- gb2fasta(plik_gbff, paste0(plik_gbff, ".fasta"));
  plik_fasta <- read.fasta(plik_fasta);
  gbff <- readLines(plik_gbff)
  b <- str_extract_all(gbff, pattern ="16s ribosomal")
  readlines_nr <- grep(pattern = "16S ribosomal", gbff[1:length(gbff)], ignore.case = T)

  lokacja_genu <- integer()
  for(i in seq(1, 10)){
    check_if_gene <- readlines_nr - i;
    check <- gbff[check_if_gene];
    check <- grep(pattern = "gene | rRNA", check, ignore.case = T);
    if(length(check) != 0){
      lokacja_genu <- i;
      break
    }
    else{next}
  }
  
  readlines_seq_nr <- readlines_nr - lokacja_genu
  r <- gbff[readlines_nr]
  t <- sub("                     /product=\"" , "", r)
  t <- sub("\"", "", t)

  for(i in seq(1, 12)){
    check_name <- gbff[i];
    check <- grep(pattern = "  ORGANISM  ", check, ignore.case = T);
    if(length(check) != 0){
      o_n <- i;
      break
    }
    else{next}
  }
  
  org_name <- gbff[o_n]
  org_name <- sub("  ORGANISM  ", "", org_name)

  d <- gbff[readlines_seq_nr]
  f <- sub("     rRNA            ", "", d)
  f <- sub("     rRNA            ", "", f)
  f <- sub("\\(", " ", f)
  f <- sub("\\)", "", f)
  f <- sub("\\.\\.", " ", f)
  f <- str_split_fixed(f, " ", 3)
  file.path <- "/home/jarek/extracting_sequences/16S_rRNA"
  
  for(ii in (seq(1, nrow(f)))){
    if(f[ii,1] == "complement"){
      start = f[ii,2]; end = f[ii,3]; a <- comp(plik_fasta[[1]][(as.integer(end)):(as.integer(start))]);
      setwd("/home/jarek/extracting_sequences/16S_rRNA");
      write.fasta(a, paste0(org_name, ", ",  t[ii]), file.out = paste0(org_name, "_",  t[[ii]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
      setwd("/home/jarek/extracting_sequences/Fragmented");
}
    else{start = f[ii,1]; end = f[ii,2]; a <- plik_fasta[[1]][(as.integer(start)):(as.integer(end))];
      setwd("/home/jarek/extracting_sequences/16S_rRNA");
      write.fasta(a, paste0(org_name, ", ",  t[ii]), file.out = paste0(org_name, "_",  t[[ii]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
      setwd("/home/jarek/extracting_sequences/Fragmented");
}
  }
}
