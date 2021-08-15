rm(list=ls())
source("TADfit/scripts/regevaluate/regevaluate.R")

filein.y<- "TADfit/TADfit_Rproject/out/y.out"
filein.p<- "TADfit/TADfit_Rproject/out/p.out"
fileout <- "path to output/tradition_regevaluate_GM12878_150000000_151000000.TIFF"
regevaluate(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, selpiter=10, fileout)


rm(list=ls())
setwd("TADfit/scripts//regevaluate/R2vsiteration.R")

filein.y<- "TADfit/TADfit_Rproject/out/y.out"
filein.p<- "TADfit/TADfit_Rproject/out/p.out"
fileout <- "path to output/tradition_R2vsiteration_GM12878_150000000_151000000.TIFF"
R2vsiteration(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, piternum=10, fileout)
