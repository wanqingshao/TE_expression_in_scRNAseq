suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))


option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path to featurecount file"),
  make_option(c("-r", "--reference"),
              type="character",
              help="Path to saf reference file"),
  make_option(c("-n", "--name"),
              type="character",
              help="name of the output data"),
  make_option(c("-p", "--thread"),
              type="numeric",
              default = 5, 
              help="number of cores default 5"),
  make_option(c("-s", "--step"),
              default = 50,
              type="numeric",
              help="Maximum number of steps for EM algorithm default 50"),
  make_option(c("-m", "--maxn"),
              default = 1,
              type="numeric",
              help="read chagne cutoff for EM algorithm default 1"))


opt <- parse_args(OptionParser(option_list=option_list))
options(stringsAsFactors = FALSE)
options(scipen=999)

suppressPackageStartupMessages(library(magrittr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(parallel, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
setDTthreads(threads = opt$thread)

#avg_frag_length <- median(nchar(head(assigned_df, 5000)$V11))


cigar_to_pos <- function(start, cigar){
  num_list <- strsplit(cigar, "[A-Z]") %>% unlist() %>% as.numeric()
  anno_list <- strsplit(cigar, "[0-9]") %>% unlist()
  anno_list <- anno_list[anno_list != ""]
  start_list <- c()
  end_list <- c()
  for(x in 1:length(anno_list)){
    if(anno_list[x] == "M"){
      new_start <- start
      new_end <- new_start + num_list[[x]]
      start_list <- c(start_list, new_start)
      end_list <- c(end_list, new_end)
      start <- new_end
    }
    if(anno_list[x] == "N"){
      start <- start + num_list[[x]]
    }
  }
  data.table(start = start_list, end = end_list)
}


get_total_feature_length <- function(redistri_needed_features, references, thread){
  width_df <- mclapply(redistri_needed_features, function(x){
    total_width <- references[references$GeneID == x] %>% makeGRangesFromDataFrame() %>% width() %>% sum()
    data.table(feature = x, total_width =total_width )
  },mc.cores = thread) %>% do.call(rbind,.)
  width_df
}

collect_features <- function(feature){
  strsplit(feature, ",", fixed = T) %>% unlist() %>%
    sort() %>% unique() %>%
    paste(., collapse  = ",")
}


reditribute_EM <- function(quant_df, redistri_summary_n, nsteps , threads, maxn, width_df){
  quant_df_original <- quant_df
  max_change = max(100, maxn)
  redistri_df_list <- list()
  for(i in 1:nsteps){
    message(opt$file," EM",i)
    df <- mclapply(unique(redistri_summary_n$features), function(x){
      features <- strsplit(x, ",", fixed = T) %>% unlist()
      number_of_reads <- redistri_summary_n[redistri_summary_n$features == x]$N
      prior <- quant_df[feature %in% features]
      if(sum(!features %in% quant_df$feature) >0){
        prior2 <- data.table(feature = features[!features %in% quant_df$feature],N = 0)
        prior <- rbind(prior, prior2)
      }
      width_df_sub <- width_df[prior$feature, ]
      if(i == 1){
        prior$prior <- prior$N / width_df_sub$total_width
        prior$prior[is.na(prior$prior)] <- 0
        prior$new_distri <- 0
        no_uniq <- which(prior$N == 0)
        with_uniq <- which(prior$N > 0)
        prior$new_distri[no_uniq] <- 1 / nrow(prior)
        prior$new_distri[with_uniq] <- (1-length(no_uniq) / nrow(prior)) * 
          prior$prior[with_uniq] / sum(prior$prior[with_uniq])
      }else{
        prior$prior <- prior$N / width_df_sub$total_width
        prior$new_distri <- 1/ sum(prior$prior) * prior$prior
      }
      #prior$reads <- x
      prior$new_distri <- prior$new_distri * number_of_reads
      prior
    }, mc.cores = threads) %>% do.call(rbind, .)
    redistributed_reads <- df[, list(N = sum(new_distri)), by = "feature"]
    redistri_df_list <- c(redistri_df_list, list(redistributed_reads$N))
    if(i >1){
      max_change <- max(redistri_df_list[[i]] - redistri_df_list[[i-1]])
    }
    message(max_change)
    quant_df <- rbind(quant_df_original, redistributed_reads)[, list(N = sum(N)), by = "feature"]
    if(max_change < maxn){
      break
    }
  }
  quant_df
}

evenly_distribute <- function(quant_df, redistri_summary_n){
  df <- redistri_summary_n[, list(feature = unlist(strsplit(features, ",", fixed = T)), 
                                  N = N / length(unlist(strsplit(features, ",", fixed = T)))),
                           by = "features"]
  redistributed_reads <- df[, list(N = sum(N)), by = "feature"]
  quant_df <- rbind(quant_df, redistributed_reads)[, list(N = sum(N)), by = "feature"]
  quant_df
}


assigned_df <- fread(opt$file, sep = "\t")
references <- fread(opt$reference)

assigned_df$feature <- gsub(".*:", "", assigned_df$V12)
selected_cols <- c("V1", "V12", "feature")


assigned_df_simp <- assigned_df[, ..selected_cols] %>% unique()
dup_names <- assigned_df_simp$V1[duplicated(assigned_df_simp$V1)]
rm(assigned_df)

unique_reads <- assigned_df_simp[!assigned_df_simp$V1 %in% dup_names]
unique_reads <- unique_reads[grep(",", unique_reads$V12, invert = T), ]
unique_reads_quant <- unique_reads[, .N, by = "feature"]

redistri_needed <- assigned_df_simp[!assigned_df_simp$V1 %in% unique_reads$V1]
redistri_needed_features <- strsplit(redistri_needed$feature, ",") %>% unlist() %>% unique()

total_width_df <- get_total_feature_length(redistri_needed_features, references, opt$thread)
total_width_df<- as.data.frame(total_width_df)
rownames(total_width_df) <- total_width_df$feature

rm(dup_names)
rm(references)
rm(unique_reads)
gc()


if(nrow(redistri_needed) >0){
  redistri_summary <- redistri_needed[, list(features = collect_features(feature)), 
                                      by = "V1"]
  redistri_summary_n <- redistri_summary[,.N, by = "features"]
  
  EM_quant_df <-reditribute_EM(unique_reads_quant, redistri_summary_n, opt$step, opt$thread, opt$maxn, total_width_df) 
  colnames(EM_quant_df)[2] <- "EM_distri" 
  even_quant_df <- evenly_distribute(unique_reads_quant, redistri_summary_n)
  colnames(even_quant_df)[2] <- "even_distri" 
  colnames(unique_reads_quant)[2] <- "uniq"
  out_df <- merge(unique_reads_quant, EM_quant_df, all = T) %>% merge(., even_quant_df)
  out_df$uniq[is.na(out_df$uniq)] <- 0
}else{
  colnames(unique_reads_quant)[2] <- "uniq"
  out_df <- unique_reads_quant
  out_df$EM_distri <- out_df$uniq
  out_df$even_distri <- out_df$uniq
}

write.table(out_df, file = opt$name, row.names = F, col.names = T, sep = "\t", quote = F)
