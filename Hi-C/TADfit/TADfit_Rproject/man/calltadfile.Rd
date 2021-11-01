\name{calltad}
\alias{calltad}
\title{
Call hierarchical TADs on replicate Hi-C data
}
\description{
Call hierarchical TADs on replicate Hi-C data. TADs are filtered based on PM test.
The function returns hierarchical TADs,their regression coefficients and pvalues based on PM test.
}
\arguments{
  \item{dirin}{path of a directory which contains replicate Hi-C contact matrix file(s)}
  \item{resolution}{bin size of the Hi-C contact matrix}
  \item{start}{start coordinate of the region of chromosome
  }
  \item{end}{end coordinate of the region of chromosome}
  \item{species}{"hg19" for human or "mm9" for mouse}
  \item{chr}{chromosome number, e.g. "chr1". chr and species are used to checke the start and end coordinates.}
  \item{winsize}{window size when calling TAD boundaries. This parameter is from TopDom.}
  \item{minisizescale}{a factor indicates the minimum size of a TAD}
  \item{maxsizescale}{a factor indicates the maximum size of a TAD}
  \item{repeatnum}{repeat numbers when solving the model}
  \item{iteration}{iteration times when solving the model}
  \item{alpha,beta}{parameters for pre-coordinate learning rate when optimizing the model}
  \item{l}{adjust the strength of l1 regularization }
}
\usage{
calltad(dirin,resolution, start, end, 
species = "hg19", chr = "chr1", winsize = 5, 
minisizescale = 0.01, maxsizescale = 0.5, 
repeatnum = 5,iteration = 5, alpha = 0.01, 
beta = 1, l = 2)
}
\examples{
\dontrun{
library(TADfit)
dirin <- "/media/disk1/liuerhu/hicda_liuerhu/github/data/expdata/GM12878-chr1"
resolution <- 25000
start <- 0
end <- 5000000
TADs <- calltad(dirin, resolution, start, end)  #Other parameters are default
}
}
