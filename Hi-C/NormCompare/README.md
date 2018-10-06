1. Source code for six H-C normalization methods. The code has been tested under Ubuntu14.04 LTS.

 1> normSCNfun.R
 Function: Implementation of SCN method
 Depends:./sparse2matrix.R (provided)
         ./gethicobjlist.R (provided)
         ./utilities/ scriptforSCN.sh (provided)
         ./SCN_sumV2.m (provided)
         R (>= 3.5.1)
         R.matlab (R package, )
         MATALB (>= 8.6.0.267246 (R2015b))
 Example:
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normSCNfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirtemp <- ./data/temp
         normSCNfun(dirin=dirin1, dirout=dirout1, dirtemp=dirtemp, species="hg19", chr="chr1", resolution=1000000)

2> normHiCNormfun.R
 Function: Implementation of HiCNorm method
 Depends: ./sparse2matrix.R (provided)
         ./gethicobjlist.R (provided)
         R (>= 3.5.1)
         BSgenome.Hsapiens.UCSC.hg19 (R package, 1.4.0)
         BSgenome.Mmusculus.UCSC.mm9 (R package, 1.4.0)
         HiTC (R package, 1.24.0)
         Mappability files corresponding to hg19 and mm9 (Placed in extdata folder of HiTC)
 Example:
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normHiCNormfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirtemp <- ./data/temp
         enzymes <- c("MboI","MboI"," MboI"," MboI"," MboI"," MboI")
         normHiCNormfun(dirin=dirin1, dirout=dirout1, dirtemp=dirtemp, species="hg19", chr="chr1", resolution=1000000, enzymes=enzymes)

