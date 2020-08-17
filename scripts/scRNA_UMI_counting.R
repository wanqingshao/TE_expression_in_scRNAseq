suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))


option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path to featurecount file, seperate with ',' if more than one "),
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

message(opt$file)

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


reditribute_EM <- function(quant_df, redistri_summary_n, nsteps , threads, maxn, width_df, cell_group){
  quant_df_original <- quant_df
  max_change = max(100, maxn)
  redistri_df_list <- list()
  for(i in 1:nsteps){
    df <- mclapply(unique(redistri_summary_n$features), function(x){
      features <- strsplit(x, ",", fixed = T) %>% unlist()
      number_of_reads <- redistri_summary_n[redistri_summary_n$features == x]$N
      prior <- quant_df[feature %in% features,]
      if(sum(!features %in% quant_df$feature) >0){
        prior2 <- data.table(feature = features[!features %in% quant_df$feature],N = 0)
        prior <- rbind(prior, prior2)
      }
      width_df_sub <- width_df[prior$feature, ]
      if(i == 1){
        prior$prior <- prior$N / width_df_sub$total_width
        #prior$prior[is.na(prior$prior)] <- 0
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



references <- fread(opt$reference)

ffiles  <- strsplit(opt$file, ",", fixed =T)
reads_output <- data.table()

for(f in ffiles){
  reads<-data.table::fread(f, na.strings=c(""),
                           select=c(1,2,3,4),header=T,fill=T,colClasses = "character" , col.names = c("RN","RG","UB","feature") )
  reads <- unique(reads)

  multi_mapped <- reads$RN[duplicated(reads$RN)]
  multi_overlap <- reads$RN[grep(",", reads$feature)]
  redistri_needed <- reads[reads$RN %in% c(multi_mapped, multi_overlap)]
  redistri_needed_features <- strsplit(redistri_needed$feature, ",") %>% unlist() %>% unique()
  total_width_df <- get_total_feature_length(redistri_needed_features, references, opt$thread)
  total_width_df<- as.data.frame(total_width_df)
  rownames(total_width_df) <- total_width_df$feature

  uniq_reads <- reads[!reads$RN %in% c(multi_mapped, multi_overlap)]
  unique_reads_quant <- unique(uniq_reads[, 2:4])[, .N, by = c("RG","feature")]
  redistri_summary <- redistri_needed[, list(features = collect_features(feature)),
                                     by = c("RG", "UB")]
  redistri_summary_n <- redistri_summary[,.N, by = c("RG", "features")]
  rm(reads)
  rm(uniq_reads)
  rm(multi_mapped)
  rm(multi_overlap)
  rm(redistri_needed)
  rm(redistri_summary)
  gc()

  redistributed_df <- lapply(unique(redistri_summary_n$RG), function(x){
    unique_reads_quant_sub <- unique_reads_quant[RG == x, 2:3]
    redistri_summary_n_sub <- redistri_summary_n[RG == x, 2:3]

    if(nrow(redistri_summary_n_sub) >0){

      EM_quant_df <- reditribute_EM(unique_reads_quant_sub, redistri_summary_n_sub, opt$step, 1, opt$maxn, total_width_df, x)
      colnames(EM_quant_df)[2] <- "EM_distri"
      even_quant_df <- evenly_distribute(unique_reads_quant_sub, redistri_summary_n_sub)
      colnames(even_quant_df)[2] <- "even_distri"
      colnames(unique_reads_quant_sub)[2] <- "uniq"
      out_df <- merge(unique_reads_quant_sub, EM_quant_df, all = T) %>% merge(., even_quant_df)
      out_df$uniq[is.na(out_df$uniq)] <- 0
      out_df$RG <- x
    }else{
      out_df <- data.table(feature = unique_reads_quant_sub$feature,
                           uniq = unique_reads_quant_sub$N,
                           EM_distri = unique_reads_quant_sub$N,
                           even_distri = unique_reads_quant_sub$N,
                           RG = x)
    }
    out_df
  }) %>% do.call(rbind, .)

  reads_output <- rbind(reads_output, redistributed_df)
}

write.table(reads_output, file = paste0(opt$name, "_count_dt.txt"),
            row.names = F, quote = F, sep = "\t")
