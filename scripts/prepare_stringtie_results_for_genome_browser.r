suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(magrittr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(parallel, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path to stringtie file"),
  make_option(c("-n", "--name"),
              type="character",
              help="name of the output data"),
  make_option(c("-p", "--thread"),
              type="numeric",
              default = 5, 
              help="number of cores"),
  make_option(c("-s", "--stranded"),
              default = FALSE,
              action="store_true",
              help="parameter to indicate whether stranded mode should be included"))


opt <- parse_args(OptionParser(option_list=option_list))
options(stringsAsFactors = FALSE)
stringtie <- fread(opt$file, skip = 2)
stringtie$gene_id <-  gsub(".*gene_id \"","", stringtie$V9) %>% gsub("\".*", "", .)
stringtie$tx_id <-  gsub(".*transcript_id \"","", stringtie$V9) %>% gsub("\".*", "", .)

stringtie_out <- mclapply(unique(stringtie$tx_id), function(x){
  gene <- subset(stringtie, tx_id == x)
  if(opt$stranded){
    strand <- as.character(gene$V7)
  }else{
    strand <- "+"
  }
  data.table(chr = gene$V1[1], 
             transcript_start = gene$V4[1],  
             transcript_stop = gene$V5[1], 
             translation_start  = min(gene$V4), 
             translation_stop = max(gene$V5), 
             strand = unique(strand),
             gene_name = gene$gene_id[1], 
             transcript_id = gene$tx_id[1],
             type = "stringtie", 
             exon_start = paste(gene$V4[-1], collapse = ","), 
             exon_stop = paste(gene$V5[-1], collapse = ","))
}, mc.cores = opt$thread) %>% do.call(rbind, .)

write.table(stringtie_out, opt$name, quote = F, row.names = F, col.names = F, sep = "\t")
