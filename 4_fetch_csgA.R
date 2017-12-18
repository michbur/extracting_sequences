library(rentrez)
library(seqinr)
library(stringr)

#### uploading gbff and fasta files to vectors

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
  
  
#### grep /gene="csgA"
  
  # b <- str_extract_all(gbff, pattern ="/gene=\"csgA\"")
  
  readlines_nr <- grep(pattern = "/gene=\"csgA\"", gbff[1:length(gbff)], ignore.case = T)
  readlines_nr
  readlines_nr_duplikat <- grep(pattern = "/gene=\"csgA\"", gbff[1:length(gbff)], ignore.case = T)
  
  # readlines_nr <- grep(pattern = "     CDS             ", gbff[31474], ignore.case = T)
  
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
  
  druga_linijka <- integer()
  
  for(iv in seq(5, 15)){
    check_name <- gbff[iv];
    check <- grep(pattern = "Bacteria; Proteobacteria; Gammaproteobacteria;", check_name, ignore.case = T);
    if(length(check) != 0){
      if(iv - 1 == o_n){
        break
      }
      else{
        druga_linijka <<- iv - 1;
        break
      }
      break
    }
    # else{next}
  }
  
################# getting right name
  
  if(length(druga_linijka) != 0){
    org_name1 <- gbff[o_n]
    org_name2 <- gbff[druga_linijka]
    org_name <- paste(org_name1, org_name2)
  }
  if(length(druga_linijka) == 0){org_name <- gbff[o_n]}
  
  o_n
  druga_linijka
  
  org_name <- sub("  ORGANISM  ", "", org_name)
  org_name <- sub("             ", " ", org_name)
  org_name <- sub("/", " ", org_name)
  print(org_name)
  
  lokacja_genu <- integer()
  check <- 0
  gen <- 0
  #### fetching cds location which is above gene    
  for(ii in seq(1,2)){
    gen <- 0
    check_if_gene <- readlines_nr[ii] - 1;
    check <- gbff[check_if_gene];
    check <- grep(pattern = "     CDS             ", check, ignore.case = F);
    if(length(check) != 0){
      lokacja_genu <- readlines_nr[ii] - 1;
      gen <- ii
      break
    }
  }
  lokacja_genu
  gen
  
  if(length(lokacja_genu) != 0){
    
    print(lokacja_genu)
    
    
#### splitting strings to single chars
    
    readlines_seq_nr <- readlines_nr[gen]
    r <- gbff[readlines_seq_nr]
    t <- sub("                     /gene=\"" , "", r)
    t <- sub("\"", "", t)
    t

    d <- gbff[lokacja_genu]
    f <- sub("     CDS             ", "", d)
    f <- sub("     CDS             ", "", f)
    f <- sub("\\(", " ", f)
    f <- sub("\\)", "", f)
    f <- sub("\\.\\.", " ", f)
    f <- str_split_fixed(f, " ", 3)
    f
    file.path <- "/home/jarek/extracting_sequences/csgA"
    
#### gene location is marked and cut from fasta with reversing complement strings (comp is getting errors for now)  
    
    # try(
      for(n in seq(1, nrow(f))){
      if(f[n,1] == "complement"){
        start = as.integer(f[n,2]); #start = as.integer(start)
        end = as.integer(f[n,3]); #end = as.integer(end)
        #a <- comp(as.character(plik_fasta[[1]][end:start]));
        a <- comp(plik_fasta[[1]][end:start])
        setwd("/home/jarek/extracting_sequences/csgA");
        # if(length(readlines_nr) != length(readlines_nr_duplikat)){
        write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  "csgA", ".fasta"), open = "a", nbchar = 60, as.string = FALSE)
        # }
        # else{
          # write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        # }
        setwd("/home/jarek/extracting_sequences/Fragmented")
      }
      if(f[n,1] == ""){next}
      if(f[n,1] != "complement"){
        start = as.integer(f[n,1]); #start = as.integer(start)
        end = as.integer(f[n,2]); #end = as.integer(end)
        a <- plik_fasta[[1]][end:start];
        setwd("/home/jarek/extracting_sequences/csgA");
        # if(length(readlines_nr) != length(readlines_nr_duplikat)){
        write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  "csgA", ".fasta"), open = "a", nbchar = 60, as.string = FALSE)
        # }
        # else{
        #   write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        # }
        # 
        # if(f[n,1] == ""){
        #   # write.fasta(a, paste0(org_name, ", ",  t[n]), file.out = paste0(org_name, "_",  t[[n]], ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
        # }
        setwd("/home/jarek/extracting_sequences/Fragmented")
      }
    }
    # )    
  }else{next}
}

