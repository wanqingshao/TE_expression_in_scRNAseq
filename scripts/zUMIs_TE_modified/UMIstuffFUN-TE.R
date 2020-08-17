fixMissingOptions <- function(config){
  if(is.null(config$barcodes$automatic)){
    if(is.null(config$barcodes$barcode_num) & is.null(config$barcodes$barcode_file)){
      config$barcodes$automatic <- TRUE
    }else{
      config$barcodes$automatic <- FALSE
    }
  }

  if(is.null(config$barcodes$BarcodeBinning)){
    config$barcodes$BarcodeBinning <- 0
  }
  if(is.null(config$barcodes$nReadsperCell)){
    config$barcodes$nReadsperCell <- 100
  }

  if(is.null(config$barcodes$demultiplex)){
    config$barcodes$demultiplex <- FALSE
  }

  if(is.null(config$counting_opts$primaryHit)){
    config$counting_opts$primaryHit <- TRUE
  }

  if(is.null(config$counting_opts$strand)){
    config$counting_opts$strand <- 0
  }

  if(is.null(config$counting_opts$velocyto)){
    config$counting_opts$velocyto <- FALSE
  }

  if(is.null(config$counting_opts$write_ham)){
    config$counting_opts$write_ham <- FALSE
  }

  return(config)
}

splitRG<-function(bccount,mem){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  if( (maxR > 2e+09 & opt$read_layout == "SE") | (maxR > 1e+09 & opt$read_layout == "PE") ){
    maxR <- ifelse(opt$read_layout == "SE",2e+09,1e+09)
  }

  print(paste(maxR,"Reads per chunk"))
  nc<-nrow(bccount)
  cs=0
  chunkID=1
  bccount[,chunkID:=0]
  for(i in 1:nc){
    cs=cs+bccount[i]$n
    if(bccount[i]$n>maxR){
      print(paste("Warning: Barcode",bccount[i]$XC,"has more reads than allowed for the memory limit!
                  Proceeding anyway..."))
    }
    if(cs>=maxR){
      chunkID=chunkID+1
      cs=bccount[i][,"n"]
    }
    bccount[i][,"chunkID"]=chunkID
  }
  return(bccount)
}


prep_samtools <- function(featfiles,bccount,cores,samtoolsexc){
  print("Extracting reads from bam file(s)...")
  nfiles=length(featfiles)
  nchunks <- length(unique(bccount$chunkID))
  all_rgfiles <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".RGgroup.",1:nchunks,".txt")


  for(i in unique(bccount$chunkID)){
    rgfile <- all_rgfiles[i]
    chunks <- bccount[chunkID==i]$XC
    write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
  }

  headerXX <- paste( c(paste0("V",1:4)) ,collapse="\t")
  write(headerXX,"freadHeader")

  headercommand <- "cat freadHeader > "
  samcommand <- paste(samtoolsexc," view -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x XS -@")
  grepcommand <- " | cut -f1,12,13,14 | grep -e 'XT:Z:' | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f "

  outfiles_ex <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".ex.",1:nchunks,".txt")
  system(paste(headercommand,outfiles_ex,collapse = "; "))

  if(length(featfiles)==1){
    cpusperchunk <- round(cores/nchunks,0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")

    system(paste(ex_cmd,"wait"))
  }else{
    cpusperchunk <- round(cores/(2*nchunks),0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")

    outfiles_in <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".in.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_in,collapse = "; "))

    in_cmd <- paste(samcommand,cpusperchunk,featfiles[2],grepcommand,all_rgfiles,">>",outfiles_in," & ",collapse = " ")

    system(paste(ex_cmd,in_cmd,"wait"))
  }
  system("rm freadHeader")
  system(paste("rm",all_rgfiles))

  return(outfiles_ex)
}
