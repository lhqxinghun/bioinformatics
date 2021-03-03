rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/lower.tri")
source("lower.tri.R")
filein.y<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/y.out"
filein.p<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/p.out"
fileout <- "/home/biology/Hic/hicda_liuerhu/lower.tri/figure/nolog_tradition_heatmap_GM12878_150000000_151000000.TIFF"


lower.tri(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, selpiter=10, fileout)
