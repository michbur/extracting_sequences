library(dplyr)
library(pbapply)
library(seqinr)
library(stringr)

genomes_path <- "/home/michal/Dokumenty/Genomy_backup/"

genome_files <- list.files(genomes_path)
gbff_file <- genome_files[1]

get_sequence <- function(ith_file, seq_name) {
  all_lines <- readLines(paste0(genomes_path, ith_file))
  
  #grep(pattern = "LOCUS       ", all_lines, ignore.case = TRUE)
  
  locus_ids <- grep("LOCUS", all_lines)
  
  locus_id <- matrix(c(locus_ids,
                       locus_ids[-1] - 1, length(all_lines)),
                     nrow = 2, byrow = TRUE)
  
  all_loci <- lapply(1L:ncol(locus_id), function(ith_locus) {
    all_lines[locus_id[1, ith_locus]:locus_id[2, ith_locus]]
  })
  
  lapply(all_loci, function(ith_loci_lines) {

    if(length(grep(pattern = seq_name, ith_loci_lines)) > 1) {
      genome_start <- grep(pattern = "^        1", ith_loci_lines)
      
      gene_pos <- ith_loci_lines[sapply(grep(pattern = seq_name, ith_loci_lines), function(i)
        grep("rRNA", ith_loci_lines[(i - 10):i]) + i - 11)] %>% 
        strsplit("            ") %>% 
        sapply(last) %>% 
        data.frame(raw = ., stringsAsFactors = FALSE) %>% 
        mutate(complement = grepl("complement", raw),
               x1 = as.numeric(sub(pattern = "complement(", replacement = "", 
                                   sapply(strsplit(raw, "..", fixed = TRUE), first),
                                   fixed = TRUE)),
               x2 = as.numeric(sub(pattern = ")", replacement = "", 
                                   sapply(strsplit(raw, "..", fixed = TRUE), last),
                                   fixed = TRUE)))
      
      genes <- lapply(1L:nrow(gene_pos), function(ith_row) {
        
        x1 <- gene_pos[ith_row, "x1"]
        x2 <- gene_pos[ith_row, "x2"]
        
        # x1 and x2 converted to ids in the cutted sequence
        xx1 <- x1 - floor((x1 - 1)/60)*60 
        xx2 <- x2 - floor((x1 - 1)/60)*60
        
        ith_loci_lines[(genome_start + floor((x1 - 1)/60)):(genome_start + floor((x2 - 1)/60))] %>% 
          strsplit("1 ") %>% 
          lapply(last) %>% 
          paste0(collapse = "") %>% 
          gsub(" ", "", .) %>% 
          substr(xx1, xx2) %>% 
          strsplit(split = "") %>% 
          unlist
      })
      
      genes[gene_pos[["complement"]]] <- lapply(genes[gene_pos[["complement"]]], function(ith_gene)
        rev(comp(ith_gene)))
      
      org_name <- ith_loci_lines[grep(pattern = "ORGANISM", ith_loci_lines)] %>% 
        strsplit("ORGANISM[ ]*") %>% 
        sapply(last)
      
      gene_names <- ith_loci_lines[sapply(grep(pattern = seq_name, ith_loci_lines), function(i)
        grep("gene=", ith_loci_lines[(i - 10):i])[1] + i - 11)] %>% 
        strsplit('\"', fixed = TRUE) %>% 
        sapply(last)
      
      locus_names <- ith_loci_lines[sapply(grep(pattern = seq_name, ith_loci_lines), function(i)
        grep("locus_tag=", ith_loci_lines[(i - 10):i])[1] + i - 11)] %>% 
        strsplit('\"', fixed = TRUE) %>% 
        sapply(last)
      
      # browser()
      # sapply(grep(pattern = seq_name, ith_loci_lines), function(i)
      #   ith_loci_lines[(i - 10):i])
      
      names(genes) <- paste0(">", org_name, "|", gene_names, "|", locus_names)
      
      sapply(grep(pattern = seq_name, ith_loci_lines), function(i)
        ith_loci_lines[(i - 10):i])
      
      genes
    }
  }) %>% 
    unlist(recursive = FALSE)
}

seq_rna <- pblapply(genome_files[1L:5], function(ith_file) 
  get_sequence(ith_file, "16S ribosomal RNA")) %>% 
  unlist(recursive = FALSE)

