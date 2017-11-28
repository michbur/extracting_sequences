library(rentrez)
library(seqinr)
library(stringr)

setwd("/home/jarek/extracting_sequences/Fragmented")
lista_plikow <- list.files("/home/jarek/extracting_sequences/Fragmented")
lista_plikow

####  converting to fasta, possible manual corrections to fix gb2fasta locus grep
for(i in lista_plikow){
  plik_gbff <- i;
  nazwa_fasta <- paste0(plik_gbff, ".fasta");
  print(plik_gbff);
  gb2fasta(plik_gbff, paste0("/home/jarek/extracting_sequences/Fasta/", nazwa_fasta));
  lista_plikow <- setdiff(lista_plikow, i)
}

lista_fasta <- list.files("/home/jarek/extracting_sequences/Fasta")
lista_plikow <- list.files("/home/jarek/extracting_sequences/Fragmented")

file_len <- length(lista_plikow)

### rest
for(i in seq(1, file_len)){
  setwd("/home/jarek/extracting_sequences/Fasta")
  plik_fasta <- lista_fasta[i];
  plik_fasta <- read.fasta(plik_fasta);
  setwd("/home/jarek/extracting_sequences/Fragmented")
  plik_gbff <- lista_plikow[i];
  gbff <- readLines(plik_gbff);
  
  b <- str_extract_all(gbff, pattern ="16s ribosomal")
  readlines_nr <- grep(pattern = "16S ribosomal", gbff[1:length(gbff)], ignore.case = T)
  readlines_nr
  
  o_n <- integer()
    
  for(iii in seq(7, 15)){
    check_name <- gbff[iii];
    check <- grep(pattern = "  ORGANISM  ", check, ignore.case = T);
    if(length(check) != 0){
      o_n <- iii;
      break
    }
    # else{next}
  }

  org_name <- gbff[o_n]
  org_name <- sub("  ORGANISM  ", "", org_name)
  
  lokacja_genu <- integer()
  
  for(ii in seq(1, 10)){
    check_if_gene <- readlines_nr - ii;
    check <- gbff[check_if_gene];
    check <- grep(pattern = "gene | rRNA", check, ignore.case = T);
    if(length(check) != 0){
      lokacja_genu <- ii;
      
      readlines_seq_nr <- readlines_nr - lokacja_genu
      r <- gbff[readlines_nr]
      t <- sub("                     /product=\"" , "", r)
      t <- sub("\"", "", t)
      
      d <- gbff[readlines_seq_nr]
      f <- sub("     rRNA            ", "", d)
      f <- sub("     rRNA            ", "", f)
      f <- sub("\\(", " ", f)
      f <- sub("\\)", "", f)
      f <- sub("\\.\\.", " ", f)
      f <- str_split_fixed(f, " ", 3)
      print(f)
      print(f[2,1])
      file.path <- "/home/jarek/extracting_sequences/16S_rRNA"
      
      for(n in (seq(1, nrow(f)))){
        if(f[n,1] == "complement"){
          start = f[n,2]; end = f[n,3]; a <- comp(plik_fasta[[1]][(as.integer(end)):(as.integer(start))]);
          setwd("/home/jarek/extracting_sequences/16S_rRNA");
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
          setwd("/home/jarek/extracting_sequences/Fragmented");
    }
        else{start = f[n,1]; end = f[n,2]; a <- plik_fasta[[1]][(as.integer(start)):(as.integer(end))];
          setwd("/home/jarek/extracting_sequences/16S_rRNA");
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE);
          setwd("/home/jarek/extracting_sequences/Fragmented");
    }
      }
    }
    else{next}
  }
}
