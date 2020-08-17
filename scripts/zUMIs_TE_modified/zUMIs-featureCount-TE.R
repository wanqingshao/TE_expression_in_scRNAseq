#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)

##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <-read_yaml(myYaml)
setwd(opt$out_dir)
#try(unixtools::set.tempdir(opt$out_dir))
source(paste0(opt$zUMIs_directory,"/runfeatureCountFUN-TE.R"))
source(paste0(opt$zUMIs_directory,"/featureCounts.R"))
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN-TE.R"))
options(datatable.fread.input.cmd.message=FALSE)
print(Sys.time())

samtoolsexc <- opt$samtools_exec
data.table::setDTthreads(threads=opt$num_threads)
fcounts_clib <- paste0(opt$zUMIs_directory,"/fcountsLib")



opt <- fixMissingOptions(opt)
#######################################################################
#######################################################################
##### Barcode handling & chunking

#read file with barcodecounts

#check if binning of adjacent barcodes should be run
if(opt$barcodes$BarcodeBinning > 0){
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}else{
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))
}
bccount<-splitRG(bccount=bccount, mem= opt$mem_limit)

##############################################################
##### featureCounts
saf<- read.table(opt$reference$SAF_file, header = T)
abamfile<-paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam")

fnex<-.runFeatureCount(abamfile,
                       saf=saf,
                       strand=opt$counting_opts$strand,
                       type="ex",
                       countMultiMappingReads = opt$counting_opts$countMultiMappingReads,
                       fracOverlap = opt$counting_opts$fracOverlap,
                       cpu = opt$num_threads,
                       mem = opt$mem_limit,
                       fcounts_clib = fcounts_clib)
ffiles<-fnex

if(is.null(opt$mem_limit)){
  mempercpu <- round(100/opt$num_threads/2,0)
}else{
  mempercpu <- round(opt$mem_limit/opt$num_threads/2,0)
  if(mempercpu==0){
    mempercpu <- 1
  }
}

system(paste0("for i in ",paste(ffiles,collapse=" ")," ; do ",samtoolsexc," sort -t BC -n -O 'BAM' -@ ",round(opt$num_threads/2,0)," -m ",mempercpu,"G -o $i $i.tmp & done ; wait"))
system(paste("rm",paste0(ffiles,".tmp",collapse=" ")))

samouts <- prep_samtools(featfiles = ffiles,
                         bccount   = bccount,
                         cores     = opt$num_threads,
                         samtoolsexc=samtoolsexc)

lapply(samouts, function(x){
  system(paste0("mv ", x, " ", opt$out_dir, "/",  opt$project, "_featurecount.sam.txt"))
})

print(Sys.time())
print(paste("zUMIs done"))
q()
