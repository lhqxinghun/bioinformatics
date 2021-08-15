rm(list=ls())
source("scripts/lower.tri/lower.tri.R")
filein.y<- "TADfit/TADfit_Rproject/out/y.out"
filein.p<- "TADfit/TADfit_Rproject/out/p.out"
fileout <- "path to output/nolog_tradition_heatmap_GM12878_150000000_151000000.TIFF"

lower.tri(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, selpiter=10, fileout)
