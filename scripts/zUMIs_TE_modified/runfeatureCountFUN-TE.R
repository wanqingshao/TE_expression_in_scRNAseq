suppressMessages(require(dplyr))
suppressWarnings(suppressMessages(require(GenomicRanges)))
suppressWarnings(suppressMessages(require(GenomicFeatures)))
suppressWarnings(suppressMessages(require(GenomicAlignments)))
suppressWarnings(suppressMessages(require(AnnotationDbi)))

.runFeatureCount<-function(abamfile,RG,saf,strand,type,countMultiMappingReads,fracOverlap, cpu,mem,fcounts_clib){
  print(paste0("Assigning reads to features (",type,")"))
     fc.stat<-featureCounts(files=abamfile,
                            annot.ext=saf,
                            isGTFAnnotationFile=F,
                            primaryOnly=!(countMultiMappingReads),
                            allowMultiOverlap = T,
                            fracOverlap= fracOverlap,
                            countMultiMappingReads=countMultiMappingReads,
                            fraction = T,
                            nthreads=cpu,
                            reportReads="BAM",
                            strandSpecific=strand,
                            isPairedEnd=T,
                            countChimericFragments=F,
                            fcounts_clib = fcounts_clib)$stat
     
  fn<-paste0(abamfile,".featureCounts.bam")
  nfn<-paste0(abamfile,".",type,".featureCounts.bam")

  system(paste0("mv ",fn," ",nfn,".tmp"))

  invisible(suppressWarnings(suppressMessages(gc(verbose=F))))
  return(nfn)
}
