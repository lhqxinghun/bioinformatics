\name{sparse2matrix}
\alias{sparse2matrix}
\title{
Convert a sparsed Hi-C contact matrix with three colomns to a square Hi-C contact matrix.
}
\description{
This function converts a sparsed Hi-C contact matrix with three colomns to a square Hi-C contact matrix.
}
\arguments{
  \item{sparmat}{sparse contact matrix in a triple format}
  \item{species}{"hg19" for human or "mm9" for mouse}
  \item{chr}{chromosome number, e.g. "chr1"}
  \item{resolution}{bin size of the contact matrix}
}
\usage{
sparse2matrix(sparmat, species, chr, resolution)	
}
