tadfit <- function(matlistcon1=NULL, matlistcon2=NULL, namecon1="Context1", namecon2="Context2", 
                   parnorm=list(bmedifilter=FALSE, bicnorm=FALSE, maxiter=50),
                   parhierTADs=list(winsize=5, bstatfilter=TRUE, bTADsfilter=TRUE, minsizescale=0.01, maxsizescale=0.5), 
                   parFTRLreg=list(repeatnum=3, iteration=2, lambda=0.1, alpha=c(0.0025,0.005,0.010), beta=c(1,1,1), l1=c(0.5,1,2), l2=c(0.5,1,2)))
{
  if(!requireNamespace("stats", quietly = TRUE))
  {
    stop("No package named stats")
    return(NULL)
  }
  
  parFTRLreg$lambda <- unique(parFTRLreg$lambda)
  parFTRLreg$alpha <- unique(parFTRLreg$alpha)
  parFTRLreg$beta <- unique(parFTRLreg$beta)
  parFTRLreg$l1 <- unique(parFTRLreg$l1)
  parFTRLreg$l2 <- unique(parFTRLreg$l2)

  output <- .Call(`_TADfit_rundiffTADs`, matlistcon1, matlistcon2, namecon1, namecon2, parnorm, parhierTADs, parFTRLreg)

  bcon1 <- grepl("B_Con1", colnames(output$TADs))
  bcon2 <- grepl("B_Con2", colnames(output$TADs))
  
# two context
  if(length(which(bcon1==TRUE)) >= 1 && length(which(bcon2==TRUE)) >= 1)
  {
    avedifffun <- function(x)
    {
      avediff <- mean(x[bcon1]) - mean(x[bcon2])
      return (avediff)
    }
    avediff <- apply(output$TADs, 1, avedifffun)
    output$TADs <- cbind(output$TADs, "avediff"=avediff)    
  }
  
# two context and repeatnum bigger than 2
  if(length(which(bcon1==TRUE)) >= 2 && length(which(bcon2==TRUE))  >= 2)
  {
    pvaluefun <- function(x)
    {
      # p_value <- t.test(x[bcon1], x[bcon2], paired = TRUE)$p.value
      p_value <-stats::wilcox.test(x[bcon1], x[bcon2], aired = TRUE, exact=TRUE)$p.value
      return (p_value)
    }
    p_value <- apply(output$TADs, 1, pvaluefun)
    # q_value <- qvalue(p_value, pi0=1)$qvalues
    # q_value <- qvalue(p_value)$qvalues
    # output$TADs <- cbind(output$TADs, "p-value"=p_value, "q-value"=q_value)
    output$TADs <- cbind(output$TADs, "p-value"=p_value)
  }
  else if(length(which(bcon1==TRUE))  >= 2)
  {
    # pvaluefun <- function(x)
    # {
    #   p_value <- t.test(x[bcon1])$p.value
    #   return (p_value)
    # }
    # p_value <- apply(output$TADs, 1, pvaluefun)
    # q_value <- qvalue(p_value, pi0=1)$qvalues
    # q_value <- qvalue(p_value)$qvalues
    # output$TADs <- cbind(output$TADs, "p-value"=p_value, "q-value"=q_value)
    # output$TADs <- cbind(output$TADs, "p-value"=p_value)
  }
  else if(length(which(bcon2==TRUE))  >= 2)
  {
    pvaluefun <- function(x)
    {
      p_value <- stats::t.test(x[bcon2])$p.value
      return (p_value)
    }
    p_value <- apply(output$TADs, 1, pvaluefun)
    # q_value <- qvalue(p_value, pi0=1)$qvalues
    # q_value <- qvalue(p_value)$qvalues
    # output$TADs <- cbind(output$TADs, "p-value"=p_value, "q-value"=q_value)
    output$TADs <- cbind(output$TADs, "p-value"=p_value) 
  }
  return(output)
}

