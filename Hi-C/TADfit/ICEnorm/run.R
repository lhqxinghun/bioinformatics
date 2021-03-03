####################### 50K 全段标准化

#resolution：50K
#dirin："/media/disk2/dump/GM12878/chr1/50K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE"
#startcor: 0
#endcor:NULL

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#全段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/50K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=50000, startcor=0, endcor=NULL)

#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE", 4, species="hg19", chr="chr1", resolution=50000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/50K/ICE-CPM")



####################### 100K 全段标准化

#resolution：100K
#dirin："/media/disk2/dump/GM12878/chr1/100K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/100K/ICE"
#startcor: 0
#endcor:NULL

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#全段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/100K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/100K/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=100000, startcor=0, endcor=NULL)

#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/100K/ICE", 4, species="hg19", chr="chr1", resolution=100000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/100K/ICE-CPM")


####################### 250K 全段标准化

#resolution：250K
#dirin："/media/disk2/dump/GM12878/chr1/250K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/250K/ICE"
#startcor: 0
#endcor:NULL

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#全段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/250K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/250K/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=250000, startcor=0, endcor=NULL)

#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/250K/ICE", 4, species="hg19", chr="chr1", resolution=250000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/250K/ICE-CPM")


####################### 500K 全段标准化

#resolution：500K
#dirin："/media/disk2/dump/GM12878/chr1/500K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/500K/ICE"
#startcor: 0
#endcor:NULL

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#全段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/500K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/500K/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=500000, startcor=0, endcor=NULL)

#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/500K/ICE", 4, species="hg19", chr="chr1", resolution=500000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/500K/ICE-CPM")


####################### 1M 全段标准化

#resolution：1M
#dirin："/media/disk2/dump/GM12878/chr1/1M"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/1M/ICE"
#startcor: 0
#endcor:NULL

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#全段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/1M", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/1M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=1000000, startcor=0, endcor=NULL)

#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/1M/ICE", 4, species="hg19", chr="chr1", resolution=1000000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/1M/ICE-CPM")

####################### 5K 分段标准化

#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-170M/ICE"
#startcor: 150000000
#endcor:170000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-170M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=150000000, endcor=170000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-170M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-170M/ICE-CPM")


#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/170M-200M/ICE"
#startcor: 170000000
#endcor:200000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/170M-200M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=170000000, endcor=200000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/170M-200M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/170M-200M/ICE-CPM")



#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/0M-50M/ICE"
#startcor: 0
#endcor:50000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/0M-50M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=0, endcor=50000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/0M-50M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/0M-50M/ICE-CPM")

#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/50M-100M/ICE"
#startcor: 50000000
#endcor:100000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/50M-100M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=50000000, endcor=100000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/50M-100M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/50M-100M/ICE-CPM")


#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/100M-150M/ICE"
#startcor: 100000000
#endcor:150000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/100M-150M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=100000000, endcor=150000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/100M-150M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/100M-150M/ICE-CPM")



#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-200M/ICE"
#startcor: 150000000
#endcor:200000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-200M/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=150000000, endcor=200000000)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-200M/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/150M-200M/ICE-CPM")



#resolution：5K
#dirin："/media/disk2/dump/GM12878/chr1/5K"
#dirout："/media/disk2/hic-TADrecognition/GM12878/chr1/5K/200M-249M/ICE"
#startcor: 200000000
#endcor: 249000000

#删除临时文件
file.remove(list.files("/media/disk2/data/temp",full.names = TRUE))

#加载标准化函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("sparse2matrix.R")
source("submatrix.R")
source("gethicobjlist.R")
source("zerofilter.R")
source("zerorecover.R")
source("normICEfun.R")

#分段标准化
hicnormalized <- normICEfun(dirin="/media/disk2/dump/GM12878/chr1/5K", dirout="/media/disk2/hic-TADrecognition/GM12878/chr1/5K/200M-end/ICE", dirtemp="/media/disk2/data/temp", species="hg19", chr="chr1", resolution=5000, startcor=200000000, endcor=NULL)


#加载CPM函数
rm(list=ls())
setwd("/home/biology/Hic/hicda_liuerhu/ICEnorm")
source("CPM.R")

CPM("/media/disk2/hic-TADrecognition/GM12878/chr1/5K/200M-end/ICE", 4, species="hg19", chr="chr1", resolution=5000, bsparse=FALSE, bcpm=TRUE, bmodified=TRUE, "/media/disk2/hic-TADrecognition/GM12878/chr1/5K/200M-end/ICE-CPM")


