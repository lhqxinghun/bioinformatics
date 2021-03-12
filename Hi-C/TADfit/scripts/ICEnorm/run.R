#resolution：50K
#dirin："/media/disk2/dump/GM12878/chr1/50K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE"
#startcor: 0
#endcor:NULL


file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/50K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=50000, startcor=0, endcor=NULL)

rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE", 4, species="hg19", chr="chr1", resolution=50000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE-CPM")
