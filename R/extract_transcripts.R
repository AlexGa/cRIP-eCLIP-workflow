

gtf_file <- snakemake@input[["gtf"]]
genome_fa <- snakemake@input[["genome_fa"]]
output <- snakemake@output[["tx_fa"]]

extractSeqsFromAnnotation <- function(annotation, genome_fasta, outputFile = NULL, format = "gtf", metaFeature = "tx", feature = "exon", concatenate = T, replaceUnstranded = "+"){

  if(length(annotation) == 1){
    if(file.exists(annotation)){
      print(paste0("Annotation seems to be a file. Annotation will be read from ",annotation,"."))
      annotation <- rtracklayer::import(annotation,format = format)
    }
  }
  if(!is.element(class(annotation)[1], "GRanges")){
    stop("It seems that reference does not belong to GRanges. Please refer to an annotation file (gtf, gff3, bed, ...) or provide a GRanges object.")
  }

  if(!file.exists(genome_fasta)){
    stop(paste0("Could not find genome FASTA file: ",genome_fasta))
  }  
  
  # To Do: Check if index already exists before creating
  if(!file.exists(paste0(genome_fasta,".fai"))){
    print(paste0("Generate index file for ", genome_fasta ))
    tryCatch(Rsamtools::indexFa(genome_fasta), error = function(e){
      print(paste0("Rsamtools was not able to generate a fasta index --> Try to generate index with samtools faidx"))
      system(paste("samtools faidx ", genome_fasta))
    })
  }
  fa_seqs <- Rsamtools::FaFile(genome_fasta)
  
  # if strand information is missing replace with "+" or "-", defined in replaceUnstranded
  if(sum(as.character(annotation@strand@values ) %in% c("+","-")) == 0){
    if(!is.element(replaceUnstranded, c("+","-"))){
      stop("The parameter \"replaceUnstranded\" can only assigned to the values \"+\" or \"-\"")
    }
    annotation@strand[which(as.character(annotation@strand) %in% "*")] <- replaceUnstranded
  }
  
  gr_db <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = annotation))
  
  trans_list <- c()
  
  if(feature == "gene"){
    trans_list <- GenomicFeatures::transcriptsBy(gr_db, by=feature, use.names = F)
  }else if(feature == "exon"){
    trans_list <- GenomicFeatures::exonsBy(gr_db, by=metaFeature, use.names = T)
  }else if(feature == "CDS"){
    trans_list <- GenomicFeatures::cdsBy(gr_db, by=metaFeature, use.names = T)
  }else if(feature == "intron"){
    trans_list <- GenomicFeatures::intronsByTranscript(gr_db, use.names = T)
  }else if(feature == "5UTR"){
    trans_list <- GenomicFeatures::fiveUTRsByTranscript(gr_db, use.names = T)
  }else if(feature == "3UTR"){
    trans_list <- GenomicFeatures::threeUTRsByTranscript(gr_db, use.names = T)
  }else{
    stop("Please use one of the following features: \"gene\", \"tx\", \"exon\", \"CDS\", \"5UTR\", \"3UTR\" or \"intron\".")
  }
  
  seqs <- Rsamtools::getSeq(x = fa_seqs, unlist(trans_list, use.names=TRUE), as.character=TRUE)
  
  tx_table <- setNames(object = GenomicRanges::width(trans_list@partitioning), nm = trans_list@partitioning@NAMES)
  
  elt <- rep(names(tx_table), tx_table) 
  
  seq_list <- c()
  
  if(concatenate){
    seq_list <- lapply(split(seqs, elt), paste, collapse="")
  }else{
    seq_list <- split(seqs, elt)
  }
  
  if(!is.null(outputFile)){
    seqinr::write.fasta(sequences = seq_list, file.out = outputFile, names = names(seq_list), nbchar = 80, as.string = T)
  }
  
  RSQLite::dbDisconnect(gr_db$conn)
  return(seq_list)
}

tx_seqs <- extractSeqsFromAnnotation(annotation = gtf_file, genome_fasta = genome_fa, outputFile = output,
  format = "gtf", metaFeature = "tx", feature = "exon",
  concatenate = T, replaceUnstranded = "*")