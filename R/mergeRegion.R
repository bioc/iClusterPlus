
mergeRegion <- function(peak.ct, rm.prob=0.75,rm.isolate=FALSE, min.gapwidth=10000){
  requireNamespace(c("Matrix","MatrixGenerics","GenomicRanges"))
  if(rm.prob < 0 | rm.prob > 1){
    stop("rm.prob must be greater than 0 and less than 1.")
  }
  rvar <- rowVars(peak.ct)
  dtMatrix <- peak.ct[rvar > quantile(rvar,probs=rm.prob),]
  #dtMatrix <- peak.ct
  #rows of dtMatrix are chromosome location with format chrx:123-456 or chrx-123-456
  chrloc <- do.call(rbind, strsplit(rownames(dtMatrix), "[:-]"))
  chrloc <- as.data.frame(chrloc, stringsAsFactors = FALSE)
  colnames(chrloc) <- c("chr", "start", "end")
  
  chrloc$start <- as.numeric(chrloc$start)
  chrloc$end   <- as.numeric(chrloc$end)
  
  gr <- makeGRangesFromDataFrame(chrloc, seqnames.field = "chr")
  
  reduced <- reduce(gr, min.gapwidth=min.gapwidth)
  df.reduced <- as.data.frame(reduced)
  red.chrloc <- paste0(df.reduced$seqnames,":",df.reduced$start,"-",df.reduced$end)
  
  # identify overlapped regions 
  ov <- findOverlaps(reduced, gr)
  
  group <- factor(queryHits(ov))
  
  # Build mapping matrix
  G <- sparse.model.matrix(~ 0 + group)   # original genomic regions × merged regions
  
  # Sum counts per merged regions
  S <- t(G) %*% dtMatrix       # merged regions × cells
  
  # Number of original regions per merged regions
  n <- colSums(G)
  
  # Merged regions-level averages
  region.m <- S / as.numeric(n)
  rownames(region.m) <- red.chrloc
  
  dtMerged <- NA
  if(rm.isolate == TRUE){
    dtMerged <- region.m[n>1,]
  }else {
    dtMerged <- region.m
  }
  dtMerged
}
