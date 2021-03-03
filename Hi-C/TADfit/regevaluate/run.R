rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/regevaluate")
source("regevaluate.R")

filein.y<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/y.out"
filein.p<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/p.out"
fileout <- "/home/biology/Hic/hicda_liuerhu/regevaluate/figure/tradition_regevaluate_GM12878_150000000_151000000.TIFF"
regevaluate(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, selpiter=10, fileout)


rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/regevaluate")
source("R2vsiteration.R")

filein.y<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/y.out"
filein.p<- "/home/biology/Hic/hicda_liuerhu/DiffTADs/out/p.out"
fileout <- "/home/biology/Hic/hicda_liuerhu/regevaluate/figure/tradition_R2vsiteration_GM12878_150000000_151000000.TIFF"
R2vsiteration(filein.y, filein.p, repnum=4, selycon=1, selyrep=3, selpcon=1, conum=1, piternum=10, fileout)
