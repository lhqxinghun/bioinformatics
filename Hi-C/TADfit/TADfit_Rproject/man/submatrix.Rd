\name{submatrix}
\alias{submatrix}
\title{Get a submatrix of a given Hi-C contact matrix
}
\description{Return a submatrix of the input one from the start to the end.
}
\arguments{
  \item{mat}{input Hi-C contact matrx}
  \item{species}{"hg19" for human or "mm9" for mouse}
  \item{chr}{chromosome number, e.g. "chr1"}
  \item{start}{start coordinate of the region of chromosome}
  \item{end}{end coordinate of the region of chromosome}
  \item{resolution}{bin size of the corresponding Hi-C contact matrix}
}
\usage{
submatrix(mat, species, chr, start, end, resolution)
}
\examples{
\dontrun{
library(TADfit)
Matriceslist <- loadmat("/media/disk1/liuerhu/hicda_liuerhu/github/data/expdata/GM12878-chr1")
Submat <- submatrix(Matriceslist$listcon1[[1]], "hg19", "chr1", 2500000, 5000000,25000)
}
}
