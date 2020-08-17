list.of.packages <- c("ggplot2", "ggrepel", "ggpubr", "optparse", "magrittr", "data.table", "parallel")
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
  make_option(c("-f", "--file"),
              type="character",
              help="Path to assembled transcript"),
  make_option(c("-g", "--gene"),
              type="character",
              help="Path to gene saf transcript"),
  make_option(c("-a", "--ncrna"),
              type="character",
              help="Path to ncRNA saf transcript"),
  make_option(c("-m", "--mito"),
              type="character",
              default = "unknown",
              help="Path to mitochondria saf transcript"),
  make_option(c("-r", "--repeats"),
              type="character",
              help="Path to repeat saf file"),
  make_option(c("-n", "--name"),
              type="character",
              help="basename of the output data"),
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
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggrepel, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts=F, quietly=T))

setDTthreads(threads = opt$thread)


ggplot_theme <-theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12),
        axis.title.x   = element_text(size=14, face="bold"),
        axis.title.y   = element_text(size=14, face="bold"),
        strip.text.x = element_text(size=10, face="bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.title =  element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
theme_set(ggplot_theme)

message("Loading samples")

assemble_df <- fread(opt$file)
gene_df <- fread(opt$gene)
te_df <- fread(opt$repeats)
nc_df <- fread(opt$ncrna)

get_stringtie_summary <- function(stringtie, pc_gr, te_gr, sample, reture_summary = T){
  colnames(stringtie)[c(1, 4,5, 7)] <- c("chr", "start", "end", "strand")
  stringtie$tx_id <- gsub(".*transcript_id \"", "", stringtie$V9) %>% gsub("\".*", "", .)
  stringtie_exon <- subset(stringtie, V3 == "exon")
  stringtie_exon_gr <- makeGRangesFromDataFrame(stringtie_exon, keep.extra.columns = T)
  stringtie_to_pc <- findOverlapPairs(stringtie_exon_gr, pc_gr,ignore.strand = T)
  stringtie_to_te <- findOverlapPairs(stringtie_exon_gr, te_gr, ignore.strand = T)
  
  stringtie_df <- data.table(tx_id = unique(stringtie_exon_gr$tx_id), type = "other")
  stringtie_df$type[stringtie_df$tx_id %in% stringtie_to_pc@first$tx_id] <- "exon"
  stringtie_df$type[stringtie_df$tx_id %in% stringtie_to_te@first$tx_id] <- "te"
  stringtie_df$type[(stringtie_df$tx_id %in% stringtie_to_te@first$tx_id) & (stringtie_df$tx_id %in% stringtie_to_pc@first$tx_id)] <- "exon_te"
  if(reture_summary){
    summary_df <- as.data.frame(table(stringtie_df$type))
    colnames(summary_df) <- c("type", "count")
    summary_df$sample <- sample
    summary_df$prop <- summary_df$count / sum(summary_df$count)
    summary_df$text <- paste0(summary_df$count, " (", round(summary_df$prop* 100, 1) , "%)")
    rownames(summary_df) <- summary_df$type
    summary_df <- summary_df[rev(c("exon", "exon_te","te",  "other")),]
    summary_df <- summary_df %>%
      mutate(lab.ypos = cumsum(prop) - 0.5*prop)
    
    summary_df
  }else{
    stringtie_df$sample <- sample
    stringtie_df
  }
}

get_gene_type <- function(type){
  exon_stats <- ifelse(length(grep("exon", type)), T, F)
  te_stats <- ifelse(length(grep("te", type)), T, F)
  if(te_stats & (!exon_stats)){
    gene_type <- "te"
  }
  if((!te_stats) & exon_stats){
    gene_type <- "exon"
  }
  if(te_stats & exon_stats){
    gene_type <- "exon_te"
  }
  if((!te_stats) & (!exon_stats)){
    gene_type <- "other"
  }
  gene_type
}

get_stringtie_gene_component_summary <- function(stringtie, saf_exon.gr, saf_te.gr, sample, reture_summary = T){
  component_df <-  get_stringtie_summary(stringtie,saf_exon.gr, saf_te.gr, sample, reture_summary = F)
  gene_ids <- gsub(".*gene_id \"", "", stringtie$V9) %>% gsub("\".*", "", .)
  tx_ids <- gsub(".*transcript_id \"", "", stringtie$V9) %>% gsub("\".*", "", .)
  
  gene_to_tx <- data.table(tx_id = tx_ids, gene_id = gene_ids) %>% unique()
  component_df <- merge(gene_to_tx,component_df)
  gene_component_df <- component_df[, list(type =get_gene_type(paste(unique(type), collapse = " "))), by = "gene_id"]
  if(reture_summary){
    summary_df <- as.data.frame(table(gene_component_df$type))
    colnames(summary_df) <- c("type", "count")
    summary_df$sample <- sample
    summary_df$prop <- summary_df$count / sum(summary_df$count)
    summary_df$text <- paste0(summary_df$count, " (", round(summary_df$prop* 100, 1) , "%)")
    rownames(summary_df) <- summary_df$type
    summary_df <- summary_df[rev(c("exon", "exon_te","te",  "other")),]
    summary_df <- summary_df %>%
      mutate(lab.ypos = cumsum(prop) - 0.5*prop)
    
    summary_df
  }else{
    gene_component_df$sample <- sample
    gene_component_df
  }
}  


saf_te.gr <- makeGRangesFromDataFrame(te_df, keep.extra.columns = T)
saf_pc.gr <- makeGRangesFromDataFrame(gene_df, keep.extra.columns = T)
nc_df.gr <- makeGRangesFromDataFrame(nc_df, keep.extra.columns = T)


message("QC analysis for assembled transcripts")

component_df <-  get_stringtie_gene_component_summary(assemble_df,saf_pc.gr, saf_te.gr, "assembled transcript")
component_df$type <- factor(component_df$typ, levels = c("exon",  "exon_te", "te","other"))

p1 <- ggplot(component_df, aes(x="", y=prop, fill=type))+ 
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)+
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values = c("#d18aad", "#8bc6d8","goldenrod1", "#D5D5D5")) +
  geom_text_repel(aes(y = lab.ypos, label = text),size = 3 ) +
  ggtitle("Type of assembled genes")




prepare_stringtie_exon <- function(stringtie, stranded = F){
  colnames(stringtie)[c(1, 4,5, 7)] <- c("chr", "start", "end", "strand")
  stringtie$gene_id <- gsub(".*gene_id \"", "", stringtie$V9) %>% gsub("\".*", "", .)
  tx_gr <- subset(stringtie, V3 == "exon") %>% makeGRangesFromDataFrame(., keep.extra.columns = T)
  if(!stranded){
    strand(tx_gr) <- "*"
  }else{
    tx_gr <- tx_gr[strand(tx_gr) %in% c("+", "-")]
  }
  seqlevels(tx_gr, pruning.mode  = "coarse") <- seqlevels(tx_gr) %>% grep("chr", ., value = T)
  tx_gr
}

component_info <- get_stringtie_gene_component_summary(assemble_df,saf_pc.gr, saf_te.gr, "merged", reture_summary = F)
te_only_genes <- subset(component_info, type == "te")$gene_id

assemble_exon_gr <- prepare_stringtie_exon(assemble_df, stranded = T)
assemble_exon_gr_te_only <- assemble_exon_gr[assemble_exon_gr$gene_id %in% te_only_genes]


nc_overlapping_genes <- findOverlapPairs(assemble_exon_gr_te_only, nc_df.gr, ignore.strand =)@first$gene_id %>% unique()
other_genes <- unique(assemble_exon_gr_te_only$gene_id)[!unique(assemble_exon_gr_te_only$gene_id) %in% nc_overlapping_genes]

nc_stats_df <- data.frame(type = c("nc", "other"), number = c(length(nc_overlapping_genes), length(other_genes)))
nc_stats_df$type <- factor(nc_stats_df$type,levels = c("other", "nc"))

p2<- ggplot(nc_stats_df , aes(x = "", y = number, fill = type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("peachpuff", "salmon")) +
  ggtitle("Overlapping TE only genes with ncRNA")

pdf(paste0(opt$name, "_transcript_qc.pdf"), width = 7.5, height = 3.5)
ggarrange(p1, p2, widths = c(1.5, 1))
dev.off()

message("Exporting assembled transcripts")


assemble_exon_gr_te_only_saf <- data.table(GeneID = assemble_exon_gr_te_only$gene_id, 
                                           Chr = as.character(seqnames(assemble_exon_gr_te_only)),
                                           Start = as.numeric(start(assemble_exon_gr_te_only)),
                                           End = as.numeric(end(assemble_exon_gr_te_only)),
                                           Strand = as.character(strand(assemble_exon_gr_te_only)),
                                           type = "te_tx")

if(opt$mito != "unknown"){
  chrM_saf <- fread(opt$mito)
  combined_saf <- rbind(gene_df, assemble_exon_gr_te_only_saf)
  combined_saf <- combined_saf[combined_saf$Chr != chrM_saf$Chr]
  combined_saf <- rbind(combined_saf, chrM_saf)
}else{
  combined_saf <- rbind(gene_df, assemble_exon_gr_te_only_saf)
}

write.table(combined_saf, file = paste0(opt$name, "_pc_exon_te_tx.saf"), row.names = F,
            quote=F, sep = "\t")
saveRDS(assemble_exon_gr_te_only, file = paste0(opt$name, "_te_tx.rds"))