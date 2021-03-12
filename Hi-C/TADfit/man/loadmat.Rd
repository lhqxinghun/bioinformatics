\name{loadmat}
\alias{loadmat}
\docType{package}
\title{Load Hi-C contact matrices to the R environment
}
\description{Load Hi-C contact matrices to the R environment. If isFile is set true, the dirin1 (and dirin2) should be a path to a file of Hi-C contact matrix, else it should be a directory that contains Hi-C contact matrices.
}
\arguments{
  \item{dirin1}{path of a directory of replicate Hi-C contact matrices or path of a Hi-C contact matrix file of the 1st context}
  \item{dirin2}{path of a directory or replicate Hi-C contact matrices or path of a Hi-C contact matrix file of the 2nd context}
  \item{isFile}{a logical variable indicating if the input is a file or a directory}
}
\usage{
calltad(dirin1, dirin2, isFile)	
}
\examples{
\dontrun{
library(TADfit)
Matriceslist <- loadmat("/media/disk1/liuerhu/hicda_liuerhu/github/TADfit/data/expdata/GM12878-chr1", isFile = FALSE)
}
}
