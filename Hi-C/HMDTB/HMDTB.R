# @author : Lv Hongqiang and Li Lin
# @brief : HMDTB.R is a software code to identify topological domains for given Hi-C contact matrix. 
# @version  1.0

# @fn HMDTB
# @param hicmat : matrix, matrixFile by hicdata
# @param window.size :vector, some numbers of bins to extend
#                            window.size = NULL, default value c(5£¬6£¬7)
#                            window.size = others, such as c(3, 4, 5£¬6) and so on
# @param statFilter : string, statFilter = F, not get Qvalu and not filter 
#                             statFilter = T, get Qvalu and filter
HMDTB<-function(hicmat, window.size, statFilter)
{
  disdata1 <- log10(hicmat+1)
  t1<-Sys.time()
  
  if (length(window.size) == 0)
  {
    window.size <- c(5,6,7)
  }
  else
  {
    window.size <- window.size
  }
  size_diag_max <- max(window.size)
  distance <- 15
  
  zeroM <- matrix(0, nrow=nrow(disdata1)+1, ncol=ncol(disdata1)+1)
  zeroM[2:nrow(zeroM), 2:ncol(zeroM)] <- disdata1
  
  my_matrix_integral <- get_integral(zeroM,distance)
  
  # Select the window size and do the arithmetic average
  haar_value_diag_add <- getscore(my_matrix_integral, 2)
  haar_value_diag_add <- c(matrix(0, nrow=1, ncol=length(haar_value_diag_add)))
  
  for (i in 1 : length(window.size))
  {
    haar_value_diag_num <- getscore(my_matrix_integral, window.size[i])
    haar_value_diag_add <- haar_value_diag_add + haar_value_diag_num
  }
  haar_value_diag <- haar_value_diag_add/length(window.size)
  
  x_peaks<- findPeaks(haar_value_diag) - 1

  n=length(x_peaks)
  pvalue<-getpvalue(n,disdata1,x_peaks,size_diag_max)
  
  # Caculate qvalue
  qvalue<-qvalue(pvalue)
  qvalue<-qvalue$qvalues
  
  x_peaks_resolution_start <- (x_peaks-1)*resolution
  start<- c(x_peaks_resolution_start)
  end<- c((x_peaks_resolution_start+resolution))
  pvalue<-c(pvalue)
  chrnum <- rep(chr,n )
  qvalue<-c(qvalue)
  pathname_hicdata <- strsplit(hicdata_path, "/", fixed=FALSE)
  name_hicdata <- pathname_hicdata[[1]][length(pathname_hicdata[[1]])]
  name_hicdata <- strsplit(name_hicdata, "[.]", fixed=FALSE)
  name_hicdata <- name_hicdata[[1]][1]
  # Output all maxima points
  if (statFilter == T)
  {
    outfile1 <- data.frame(chrnum, start, end, pvalue, qvalue)
    name_result <- paste(name_hicdata, "_allpeaks_byHMDTB.res", sep="")
    pathname_result <- paste("./output/", name_result,sep = "")
    write.table(outfile1, file = pathname_result, append = FALSE,
                sep = "  ", row.names = TRUE, col.names = TRUE)
    
    j=1
    boundary<-c()
    for(i in 1:n)
    {
      if(qvalue[i]< 0.05)
      {
        disdata1[x_peaks[i]+1,x_peaks[i]+1] <- 0
        boundary[j]<-x_peaks[i]
        j <- j + 1
      }
    }
    # pheatmap(disdata1,col=colorpanel(128,"lightyellow","red"),border_color=NA,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,fontsize=4)
    
    boundary_resolution <- (boundary-1)*resolution
    from.coord <- c(0, boundary_resolution)
    to.coord   <- c(boundary_resolution, (nrow(my_matrix_integral)-1)*resolution)
    start <- c(boundary_resolution)
    end <- c(boundary_resolution+resolution)
    outfile2 <- data.frame(start, end)
    name_result <- paste(name_hicdata, "_domains_byHMDTB.res", sep="")
    pathname_result <- paste("./output/", name_result,sep = "")
    
    # Output boundary point
    write.table(outfile2, file = pathname_result, append = FALSE,
                sep = "  ", row.names = TRUE, col.names = TRUE)
    list_result <- list("allpeaks" = outfile1, "domains" = outfile2)
    return(list_result)
  }
  if (statFilter == F)
  {
    outfile3 <- data.frame(chrnum, start, end)
    name_result <- paste(name_hicdata, "_allpeaks_byHMDTB.res", sep="")
    pathname_result <- paste("./output/", name_result,sep = "")
    write.table(outfile3, file = pathname_result, append = FALSE,
                sep = "  ", row.names = TRUE, col.names = TRUE)
    list_result <- list("allpeaks" = outfile3)
    return(list_result)
  }
  t2<-Sys.time()
  t=t2-t1
}

# integral graph
get_integral <- function(data,distance)
{
  
  zeroM<-data
  my_matrix_integral <- matrix(data=0,nrow=nrow(zeroM),ncol=ncol(zeroM))
  for (j in 2:(nrow(my_matrix_integral)-distance+1)) {
    for (i in j:(j+distance-1)) {
      if(i == j)
      {
        my_matrix_integral[i,j] <- 2*my_matrix_integral[i,j-1] - my_matrix_integral[i-1,j-1] + zeroM[i,j]
      }
      else if(i == (j+distance-1))
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] + zeroM[i,j]
      }
      else
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] - my_matrix_integral[i-1,j-1] + my_matrix_integral[i,j-1] + zeroM[i,j]
      }
      my_matrix_integral[j,i] <- my_matrix_integral[i,j]
    }
  }
  
  for (j in (nrow(my_matrix_integral)-distance+2):nrow(my_matrix_integral)) {
    for (i in j:nrow(my_matrix_integral)) {
      if(i == j)
      {
        my_matrix_integral[i,j] <- 2*my_matrix_integral[i,j-1] - my_matrix_integral[i-1,j-1] + zeroM[i,j]
      }
      else
      {
        my_matrix_integral[i,j] <- my_matrix_integral[i-1,j] - my_matrix_integral[i-1,j-1] + my_matrix_integral[i,j-1] + zeroM[i,j]
      }
      my_matrix_integral[j,i] <- my_matrix_integral[i,j]
    }
  }
  return(my_matrix_integral)
}

# score
getscore <- function(my_matrix_integral,size_diag)
{
  step <- 1
  i <- 1+size_diag
  j <- 1
  haar_value_diag <- matrix(data=NA,nrow=1)
  while(i <= (nrow(my_matrix_integral)-size_diag))
  {
    a1 <- my_matrix_integral[i,i]
    b1 <- my_matrix_integral[i-size_diag,i-size_diag]
    c1 <- my_matrix_integral[i,i-size_diag]
    d1 <- my_matrix_integral[i-size_diag,i]
    a2 <- my_matrix_integral[i+size_diag,i+size_diag]
    c2 <- my_matrix_integral[i+size_diag,i]
    d2 <- my_matrix_integral[i,i+size_diag]
    c3 <- my_matrix_integral[i+size_diag,i-size_diag]
    d4 <- my_matrix_integral[i-size_diag,i+size_diag]
    value <- (4*a1+b1-2*c1-2*d1+a2-2*c2-2*d2+c3+d4)/(size_diag*size_diag)
    haar_value_diag[j] <- c(value)
    j <- j+1
    i <- i+step
  }
  
  
  
  fill_diag <- matrix(0,ncol = size_diag)
  
  haar_value_diag <- c(fill_diag,haar_value_diag,fill_diag)
  haar_value_diag <- haar_value_diag/max(haar_value_diag)
  
  return(haar_value_diag)  
}

# Pvalue
getpvalue <- function(n,data,peak,size_diag)
{
  pvalue <- rep(1,n )
  left_on<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  left_down<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  right_on<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  right_down<-matrix(nrow=size_diag+1,ncol=size_diag+1)
  
  n_bins<-nrow(data)
  x_peaks<-peak
  scale.matrix.data<-hicmat
  # The data of the sub-diagonal are normalized
  for( i in 1:(2*size_diag)) 
  {
    
    scale.matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] = scale(hicmat[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ]) 
    scale.matrix.data[ seq((1+i), n_bins*(n_bins-i), 1+n_bins) ] = scale(hicmat[ seq((1+i), n_bins*(n_bins-i), 1+n_bins)]) 
  }
  
  scale.data<-matrix(NA,nrow=(n_bins+2*size_diag),ncol=(n_bins+2*size_diag))
  scale.data[(size_diag+1):(n_bins+size_diag),(size_diag+1):(n_bins+size_diag)]<-scale.matrix.data
  
  # Do rank sum test
  for(i in 1:n)
  {
    left=x_peaks[i]-size_diag
    mid=x_peaks[i]
    right=x_peaks[i]+size_diag-1
    left_on<-scale.data[left:(mid-1),left:(mid-1)]
    right_down<-scale.data[mid:right,mid:right]
    right_on<-scale.data[left:(mid-1),mid:right]
    left_down<-scale.data[mid:right,left:(mid-1)]
    left_on<-as.vector(left_on)
    left_down<-as.vector(left_down)
    right_on<-as.vector(right_on)
    right_down<-as.vector(right_down)
    wil.test =  wilcox.test(x=c(right_on,left_down), y=c(left_on,right_down), alternative="less", exact=F)
    pvalue[i] = wil.test$p.value
  }
  
  return(pvalue)
  
}

# Qvalue
qvalue <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...)
{
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  
  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }
  
  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }
  
  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
pi0est <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...) 
{
  # Check input arguments
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda) # guard against user input
  
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", ll, ".", sep=""),
                       "If length of lambda greater than 1,",
                       "you need at least 4 values.")))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  
  if (max(p) < max(lambda)) {
    stop("ERROR: maximum p-value is smaller than lambda range. Change the range of lambda or use qvalue_truncp() for truncated p-values.") 
  }
  
  # Determines pi0
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec=lambda))[ind]) / (length(p) * (1-lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      # Bootstrap method closed form solution by David Robinson
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }
  if (pi0 <= 0) {
    warning("The estimated pi0 <= 0. Setting the pi0 estimate to be 1. Check that you have valid p-values or use a different range of lambda.")
    pi0 <- pi0.lambda <- 1
    pi0Smooth <- lambda <- 0
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda,
              lambda = lambda, pi0.smooth = pi0Smooth))
}
lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...)
{
  # Check inputs
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  # Local FDR method for both probit and logit transformations
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x) / y
  } else {
    x <- log((p + eps) / (1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1 + exp(x)) ^ 2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}