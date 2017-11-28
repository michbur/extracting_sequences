library(rentrez)
library(seqinr)

setwd("/home/jarek/extracting_sequences/Fragmented")
lista_plikow <- list.files("/home/jarek/extracting_sequences/Fragmented")
lista_plikow

####  converting to fasta
for(i in lista_plikow){
  plik_gbff <- i;
  nazwa_fasta <- paste0(plik_gbff, ".fasta");
  print(plik_gbff);
  gb2fasta(plik_gbff, paste0("/home/jarek/extracting_sequences/Fasta/", nazwa_fasta));
}

lista_fasta <- list.files("/home/jarek/extracting_sequences/Fasta")
lista_fasta


### rest
for(i in length(lista_plikow)){
  plik_gbff <- lista_plikow[i];
  plik_fasta <- lista_fasta[i];
  plik_fasta <- read.fasta(plik_fasta);
  gbff <- readLines(plik_gbff)
  b <- str_extract_all(gbff, pattern ="16s ribosomal")
  readlines_nr <- grep(pattern = "16S ribosomal", gbff[1:length(gbff)], ignore.case = T)

  lokacja_genu <- integer()
  for(ii in seq(1, 10)){
    check_if_gene <- readlines_nr - ii;
    check <- gbff[check_if_gene];
    check <- grep(pattern = "gene | rRNA", check, ignore.case = T);
    if(length(check) != 0){
      lokacja_genu <- ii;
      break
    }
    else{next}
  }

  readlines_seq_nr <- readlines_nr - lokacja_genu
  r <- gbff[readlines_nr]
  t <- sub("                     /product=\"" , "", r)
  t <- sub("\"", "", t)

  for(iii in seq(7, 15)){
    check_name <- gbff[iii];
    check <- grep(pattern = "  ORGANISM  ", check, ignore.case = T);
    if(length(check) != 0){
      o_n <- iii;
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

  for(iiii in (seq(1, nrow(f)))){
    if(f[iiii,1] == "complement"){
      start = f[iiii,2]; end = f[iiii,3]; a <- comp(plik_fasta[[1]][(as.integer(end)):(as.integer(start))]);
      setwd("/home/jarek/extracting_sequences/16S_rRNA");
      write.fasta(a, paste0(org_name, ", ",  t[iiii]), file.out = paste0(org_name, "_",  t[[iiii]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
      setwd("/home/jarek/extracting_sequences/Fragmented");
}
    else{start = f[iiii,1]; end = f[iiii,2]; a <- plik_fasta[[1]][(as.integer(start)):(as.integer(end))];
      setwd("/home/jarek/extracting_sequences/16S_rRNA");
      write.fasta(a, paste0(org_name, ", ",  t[iiii]), file.out = paste0(org_name, "_",  t[[iiii]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
      setwd("/home/jarek/extracting_sequences/Fragmented");
}
  }
}
