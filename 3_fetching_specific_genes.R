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



### rest

#### uploading gbff and fasta files

lista_fasta <- list.files("/home/jarek/extracting_sequences/Fasta")
lista_plikow <- list.files("/home/jarek/extracting_sequences/Fragmented")

file_len <- length(lista_plikow)


for(i in seq(1, file_len)){
  setwd("/home/jarek/extracting_sequences/Fasta")
  plik_fasta <- lista_fasta[i];
  print(plik_fasta);
  plik_fasta <- read.fasta(plik_fasta);
  setwd("/home/jarek/extracting_sequences/Fragmented")
  plik_gbff <- lista_plikow[i];
  print(plik_gbff)
  gbff <- readLines(plik_gbff);

  
#### grep product 16S ribosomal RNA location
    
  b <- str_extract_all(gbff, pattern ="16s ribosomal RNA")
  readlines_nr <- grep(pattern = "16S ribosomal RNA", gbff[1:length(gbff)], ignore.case = T)
  readlines_nr
  readlines_nr_duplikat <- grep(pattern = "16S ribosomal RNA ", gbff[1:length(gbff)], ignore.case = F)
  #if (is.null(readlines_nr) == T ){break}
  

#### fetching organism name
  
  o_n <- integer()
  for(iii in seq(5, 15)){
    check_name <- gbff[iii];
    check <- grep(pattern = "  ORGANISM  ", check_name, ignore.case = T);
    if(length(check) != 0){
      o_n <<- iii;
      break
    }
    # else{next}
  }
  
  
  org_name <- gbff[o_n]
  org_name <- sub("  ORGANISM  ", "", org_name)
  print(org_name)
  
  lokacja_genu <- integer()
  check <- 0

  
  
#### fetching gene location which is above product location    
  for(ii in seq(1, 10)){
    check_if_gene <- readlines_nr - ii;
    check <- gbff[check_if_gene];
    check <- grep(pattern = " rRNA ", check, ignore.case = F);
    if(length(check) != 0){
      lokacja_genu <- ii;
      break
    }
  }
  if(length(lokacja_genu) != 0){
    
    # print(lokacja_genu)

    
#### splitting strings to single chars
    
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
    file.path <- "/home/jarek/extracting_sequences/16S_rRNA"

#### gene location is marked and cut from fasta with reversing complement strings (comp is getting errors for now)  
        
    for(n in seq(1, nrow(f))){
      if(f[n,1] == "complement"){print(n);print("1")
        start = f[n,2]; start = as.integer(start);print("2")
        end = f[n,3]; end = as.integer(end);print("3")
        a <- comp(as.character(plik_fasta[[1]][end:start]));print("4")
        setwd("/home/jarek/extracting_sequences/16S_rRNA");
        if(length(readlines_nr) != length(readlines_nr_duplikat)){print("nie1")
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], "_", n, ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        }
        else{print("nie2")
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        }
        setwd("/home/jarek/extracting_sequences/Fragmented")
      }
      if(f[n,1] != "complement"){print(n)
        start = f[n,1]; start = as.integer(start);
        end = f[n,2]; end = as.integer(end);
        a <- plik_fasta[start:end];
        setwd("/home/jarek/extracting_sequences/16S_rRNA");
        if(length(readlines_nr) != length(readlines_nr_duplikat)){print("tak1")
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], "_", n, ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        }
        else{print("tak2")
          write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        }
        setwd("/home/jarek/extracting_sequences/Fragmented")
      }
    }
  }else{next}
}

