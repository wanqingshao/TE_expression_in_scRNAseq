library(GenomicRanges)
library(rtracklayer)

apply_seqlengths <- function(gr, genome=Dmelanogaster) {
  seqlengths(gr) <- seqlengths(genome)[seqlevels(gr)]
  gr
}

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    if(is.null(regions)) {
      cvg <- import(cvg, as="RleList")
    } else {
      stopifnot(file.exists(cvg))
      cvg <- import(cvg, which=regions, as="RleList")
    }
  }
  cvg
}

regionApply <- function(regions, cvg, func, ...) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(unlist(
    viewApply(
      Views(cvg, as(regions, "IntegerRangesList")),
      function(x) { func(as.numeric(x), ...) },
      simplify=FALSE
    ), use.names=FALSE), use.names=FALSE)
  ans
}

regionApply_list <- function(regions, cvg, func) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewApply(
      Views(cvg, as(regions, "IntegerRangesList")),
      function(x) { func(as.numeric(x)) },
      simplify=FALSE
    ), use.names=FALSE)
  ans
}



regionSums <- function(regions, cvg, stranded = F) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewSums(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMeans <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewMeans(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewWhichMaxs(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewWhichMins(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewMaxs(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewMins(
      Views(cvg, as(regions, "IntegerRangesList"))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

total_signal <- function(cov) {
  if(class(cov) == "character") {
    cache.file <- paste0(cov, ".ts.rds")
    if(file.exists(cache.file)) {
      return(readRDS(cache.file))
    } else {
      cov <- check_coverage_argument(cov)
      ts <- sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
      saveRDS(ts, file=cache.file)
      return(ts)
    }
  } else {
    return(sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x))))))
  }
}
