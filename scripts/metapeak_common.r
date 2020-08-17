### script containing functions essential for metapeak/single-gene plot
library(GenomicRanges)
library(reshape2)
library(lattice)

read_matrix <- function(gr, cov, reverse_reads=FALSE, df=F, nu=25) {
  if(class(cov) == "character") {
    cov <- import.bw(cov, which=gr, as="RleList")
  }
  transform_function <- if(reverse_reads) { rev } else { identity }
  o <- order(gr)
  gr <- gr[o]
  
  if(df==T){
    reads.list <- regionApply_list(gr,cov, as.numeric )
    reads.list <-lapply(reads.list, function(reads) { approx(reads, n=nu)$y })
    reads.m <- matrix(as.numeric(unlist(reads.list, use.name=F)), nrow=length(reads.list), byrow=T)
    if(reverse_reads == T)reads.m <- reads.m[, ncol(reads.m):1]
  }else{
    reads.list <- regionApply(gr,cov, as.numeric )
    rl <- as(gr, "IntegerRangesList")
    view <- RleViewsList(rleList=cov[names(rl)], rangesList=rl)
    reads.list <- viewApply(view, function(x) { transform_function(as.numeric(x)) })
    reads.m <- matrix(unlist(sapply(reads.list, as.numeric)), nrow=length(gr), byrow=TRUE)
  }
  reads.m[o, ] <- reads.m
  reads.m
}

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100, diff_length = F, approx_nu=25) {
  if(diff_length == F){
    regions.gr <- resize(regions.gr, width=downstream)
    regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")
    
    reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  }else{
    reads <- matrix(nrow=length(regions.gr), ncol=approx_nu)
  }
  i_p <- which(as.vector(strand(regions.gr) == "+" | strand(regions.gr) == "*"))
  i_n <- which(as.vector(strand(regions.gr) == "-"))
  
  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")
  
  if(class(sample.cov) == "character") {
    sample.cov <- import.bw(sample.cov, which=regions.gr, as =  "RleList")
  }
  if(diff_length == F){
    if(length(i_p) > 0) reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov)
    if(length(i_n) > 0) reads[i_n, ] <- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE)
  }else{
    if(length(i_p) > 0)reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov, df=diff_length, nu=approx_nu)
    if(length(i_n) > 0)reads[i_n, ]<- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE, df=diff_length,  nu=approx_nu)
  }
  reads
}

exo_metapeak_matrix <- function(regions.gr, sample, upstream=100, downstream=100) {
  regions.gr <- resize(regions.gr, width=downstream)
  regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")
  regions.gr <- regions.gr[width(regions.gr) == upstream+downstream]
  
  i_p <- which(as.vector(strand(regions.gr) == "+" | strand(regions.gr) == "*"))
  i_n <- which(as.vector(strand(regions.gr) == "-"))
  
  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")
  
  reads.p <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  reads.n <- reads.p
  
  if(length(i_p) > 0) {
    reads.p[i_p, ] <- read_matrix(regions.gr[i_p], sample$pos, df=F, nu=NULL)
    reads.n[i_p, ] <- abs(read_matrix(regions.gr[i_p], sample$neg, df=F, nu=NULL))
  }
  
  if(length(i_n) > 0) {
    reads.p[i_n, ] <- abs(read_matrix(regions.gr[i_n], sample$neg, reverse_reads=TRUE, df=F, nu=NULL))
    reads.n[i_n, ] <- read_matrix(regions.gr[i_n], sample$pos, reverse_reads=TRUE, df=F, nu=NULL)
  }
  
  list(pos=reads.p, neg=reads.n)
}

standard_metapeak <- function(gr, sample, upstream=100, downstream=100, sample_name=NA, smooth=NA, different_length=F,approx_n=25) {
  message("standard metapeak: ", sample_name)
  if(different_length==F){
    reads <- standard_metapeak_matrix(gr, sample, upstream, downstream, diff_length=F, approx_nu=NULL)
    reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                           reads=colMeans(reads),
                           sample_name=sample_name)
    if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
    
  }else{
    reads <- standard_metapeak_matrix(regions.gr = gr, sample,upstream, downstream,  diff_length= T, approx_n=approx_n)
    reads.df <- data.frame(tss_distance=1:ncol(reads),
                           reads=colMeans(reads),
                           sample_name=sample_name)
  }
  
  reads.df
}


exo_metapeak <- function(gr, sample, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("exo metapeak: ", sample_name)
  reads.list <- exo_metapeak_matrix(gr, sample, upstream, downstream)
  
  reads.p <- reads.list$pos
  reads.n <- reads.list$neg
  
  df.p <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.p),
                     strand="+")
  
  df.n <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.n) *(-1),
                     strand="-")
  
  if(!is.na(smooth)) {
    df.n$reads <- as.numeric(runmean(Rle(df.n$reads), k=smooth, endrule="constant"))
    df.p$reads <- as.numeric(runmean(Rle(df.p$reads), k=smooth, endrule="constant"))
  }
  
  reads.df <- rbind(df.p, df.n)
  reads.df$sample_name <- sample_name
  reads.df$sample <- paste(reads.df$sample_name, reads.df$strand)
  reads.df
}
