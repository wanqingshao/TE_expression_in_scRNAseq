list.of.packages <- c("optparse", "magrittr", "data.table", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.packages2 <- c("GenomicRanges")
new.packages2 <- list.of.packages2[!(list.of.packages2 %in% installed.packages()[,"Package"])]
if(length(new.packages2)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(new.packages2)
}

suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))


option_list <- list(
  make_option(c("-g", "--gene"),
              type="character",
              help="Path to refseq annotation file"),
  make_option(c("-r", "--repeats"),
              type="character",
              help="Path to rmsk.txt file"),
  make_option(c("-n", "--name"),
              type="character",
              help="basename of the output data"),
  make_option(c("-m", "--mito"),
              type="character",
              default = "unknown",
              help="pattern for mitochondira chromosome"),
  make_option(c("-l", "--length"),
              default = "unknown",
              type="character",
              help="Path to chromosome length files"),  
  make_option(c("-p", "--thread"),
              type="numeric",
              default = 5, 
              help="number of cores default 5"))



opt <- parse_args(OptionParser(option_list=option_list))
options(stringsAsFactors = FALSE)
options(scipen=999)

suppressPackageStartupMessages(library(magrittr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(parallel, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
setDTthreads(threads = opt$thread)

message("Loading samples")

gene_df <- fread(opt$gene)
te_df <- fread(opt$repeats)
if(opt$length != "unknown" & opt$mito != "unknown"){
  chr_df <- fread(opt$length)
}


message("Generating mRNA and ncRNA saf files")
gene_df_pc <- gene_df[grep("NM_", gene_df$V2)]
gene_df_nc <- gene_df[grep("NR_", gene_df$V2)]

pc_exons <- gene_df_pc[, list(start = as.numeric(unlist(strsplit(V10, ",", fixed = T))),
                              end = as.numeric(unlist(strsplit(V11, ",", fixed = T)))),
                       by = c("V13", "V4", "V3")]
nc_exons <- gene_df_nc[, list(start = as.numeric(unlist(strsplit(V10, ",", fixed = T))),
                              end = as.numeric(unlist(strsplit(V11, ",", fixed = T)))),
                       by = c("V13", "V4", "V3")]

colnames(pc_exons) <- c("GeneID", "Strand", "Chr", "Start", "End")
colnames(nc_exons) <- c("GeneID", "Strand", "Chr", "Start", "End")

pc_exons$type <- "pc_exons"
nc_exons$type <- "nc_exons"

selected_cols <- c("GeneID", "Chr", "Start", "End","Strand",  "type")
pc_exons_saf <- pc_exons[, ..selected_cols]
nc_exons_saf <- nc_exons[, ..selected_cols]

pc_exons.gr <- makeGRangesFromDataFrame(pc_exons, keep.extra.columns = T)

message("Generating TE saf files")

te.gr <- with(te_df, GRanges(ranges = IRanges(start = V7, end = V8),
                             strand = V10, seqnames = V6,
                             repName = V11, repFamily = V13, repClass = V12))
te.gr$name <- paste0("te_", 1:length(te.gr))
te.gr <- te.gr[te.gr$repClass %in% c("DNA", "DNA?", "LINE", "LINE?", "LTR", "LTR?", "SINE", "SINE?")]


te_to_exon <- findOverlapPairs(te.gr, pc_exons.gr, ignore.strand = T)

tes_to_rm <- te_to_exon@first$name %>% unique()
te.gr_non_exon <- te.gr[!te.gr$name %in% tes_to_rm]

te_off_exon_saf <- mclapply(tes_to_rm, function(x){
  te.gr_sub <- te_to_exon@first[te_to_exon@first$name == x]
  pc_sub <- te_to_exon@second[te_to_exon@first$name == x]
  left <- setdiff(te.gr_sub, pc_sub, ignore.strand = T)
  if(length(left) >0){
    data.table(GeneID = te.gr_sub$name[1], Chr = as.character(seqnames(left)),
               Start = as.numeric(start(left)), End = as.numeric(end(left)),
               Strand = as.character(strand(te.gr_sub))[1],
               type = "te_exon_ext")
  }else{
    data.table()
  }
}, mc.cores = opt$thread) %>%do.call(rbind, .)


te_saf <- data.table(GeneID = te.gr_non_exon$name, 
                     Chr = as.character(seqnames(te.gr_non_exon)),
                     Start = start(te.gr_non_exon),
                     End = end(te.gr_non_exon),
                     Strand = as.character(strand(te.gr_non_exon)),
                     type = "te_non_exon")

te_saf <- rbind(te_saf, te_off_exon_saf)

te.gr_sub <- te.gr[te.gr$name %in% te_saf$GeneID]


if(opt$length != "unknown" & opt$mito != "unknown"){
  message("Generating chrM saf")
  chr_df <- fread(opt$length)
  chrM <- chr_df[grep(opt$mito, chr_df$V1)]
  if(nrow(chrM) >0){
    chrM_saf <- data.table(GeneID = chrM$V1, Chr = chrM$V1, Start = 1, End = chrM$V2, Strand = "+", type = "chrM")
    write.table(chrM_saf, file = paste0(opt$name, "_chrM.saf"), row.names = F, quote = F, sep = "\t")
  }else{
    message("Could not find the provided mitochondria pattern, please check the input")
  }
}

message("Saving all saf files")

write.table(te_saf, file = paste0(opt$name, "_te.saf"), row.names = F, quote = F, sep = "\t")
saveRDS(te.gr_sub, file = paste0(opt$name, "_te.gr.rds"))
write.table(pc_exons_saf, file = paste0(opt$name, "_pc_exon.saf"), row.names = F, quote = F, sep = "\t")
write.table(nc_exons_saf, file = paste0(opt$name, "_nc_exon.saf"), row.names = F, quote = F, sep = "\t")
