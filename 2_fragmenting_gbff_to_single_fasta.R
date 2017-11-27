getwd()
setwd("/home/jarek/extracting_sequences/Genomes")
lista_plików <- list.files("/home/jarek/extracting_sequences/Genomes")
lista_plików

for(i in lista_plików){
  plik_gbff <- i
  print(plik_gbff)
  multiple_enties <- readLines(plik_gbff)
  multiple_enties <- grep(pattern = "LOCUS       ", multiple_enties[1:length(multiple_enties)], ignore.case = T)
  f <- file(plik_gbff, open="rb");
  nlines <- 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
    nlines <- nlines + sum(chunk == as.raw(10L))
  }
  close(f);
  multiple_enties <- c(multiple_enties, (nlines+1))

  if(length(multiple_enties) > 2){
    for(iii in seq(1, length(multiple_enties), by = 1)){
      val <- iii+1;
      first <- multiple_enties[iii]
      qwerty <- readLines(plik_gbff);
      ostatni <- multiple_enties[val];
      
      if (val > length(multiple_enties)){break}
      else if(val <= length(multiple_enties)){
        kolejna_sekwencja <- qwerty[first:(ostatni - 1)];
        setwd("/home/jarek/extracting_sequences/Fragmented")
        fileConn <- file(paste0(plik_gbff, "_", iii, ".gbff"));
        writeLines(kolejna_sekwencja, fileConn);
        close(fileConn)
        setwd("/home/jarek/extracting_sequences/Genomes")
        print(c(i, iii))
      }
      else if(val == length(multiple_enties)){
        kolejna_sekwencja <- qwerty[first:(ostatni)];
        setwd("/home/jarek/extracting_sequences/Fragmented")
        fileConn <- file(paste0(plik_gbff, "_", iii, ".gbff"));
        writeLines(kolejna_sekwencja, fileConn);
        close(fileConn);
        setwd("/home/jarek/extracting_sequences/Genomes");
        print(c(i, iii));
        lista_plików <- setdiff(lista_plików, i)
        ### file.remove(plik_gbff)
      }
    }
  }
  else{file.copy(plik_gbff, "/home/jarek/extracting_sequences/Fragmented");
    lista_plików <- setdiff(lista_plików, i)}  
}
