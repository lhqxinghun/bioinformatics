calltad <- function(dirin, resolution, start, end, species = "hg19",  chr = "chr1", winsize = 5, minisizescale = 0.01, maxsizescale = 0.5, repeatnum = 5, iteration = 5, alpha = 0.01, beta = 1, l = 2)
{
  CPMnormatlist_read <- loadmat(dirin1=dirin, dirin2=NULL, isFile = F)
  CPMnormatlist <- lapply(CPMnormatlist_read$listcon1, submatrix, species=species, chr = chr,start, end, resolution=resolution)
  parhierTADs=list(winsize=winsize, bstatfilter=TRUE, bTADsfilter=TRUE, minsizescale=minisizescale, maxsizescale=maxsizescale)
  parFTRLreg=list(repeatnum=repeatnum, iteration=iteration, lambda=0.1, alpha=alpha, beta=beta, l1=l, l2=0)
  
  output <- tadfit(matlistcon1=CPMnormatlist, matlistcon2=NULL, parhierTADs = parhierTADs, parFTRLreg=parFTRLreg)
  
  if(requireNamespace("perm", quietly = TRUE))
  {
    pvaluefunpermTS_1_zero <- function(x)
    {
      p_value <- perm::permTS(x, 0, alternative = "greater", method = "pclt")$p.value
      return (p_value)
    }
  }
  else
  {
    stop("No package named matrixStats")
    return(NULL)
  }
  if(requireNamespace("matrixStats", quietly = TRUE))
  {
    threshold <- 0.05
    TADs <- output$TADs
    colmax <- matrixStats::colMaxs(TADs[, 5:ncol(TADs)])
    colmin <- matrixStats::colMins(TADs[, 5:ncol(TADs)])
    b_temp <- rowMeans(TADs[, 5:ncol(TADs)])
    TADs[, 5:ncol(TADs)] <- t((t(TADs[, 5:ncol(TADs)]) - colmin)/(colmax - colmin))
    blarge <- rep(TRUE, nrow(TADs))
    for (i in 5:ncol(TADs)) {
      blarge <- blarge & (TADs[,i] > threshold) 
    }
    TADs <- TADs[blarge, ]
    pvaluepm_1_zero <- apply(TADs[, 5: ncol(TADs)], 1, pvaluefunpermTS_1_zero)
    TADs <- TADs[, 1:2]
    TADs <- cbind(TADs, "b_avg" = b_temp[blarge])
    TADs <- cbind(TADs, "pvalue" = pvaluepm_1_zero)
    TADs <- TADs[TADs[, "pvalue"] < 0.05, ]
    colnames(TADs) <- c("start", "end", "b_avg", "pvalue")
  }
  else
  {
    stop("No package named perm")
    return(NULL)
  }
  return(TADs)
}
