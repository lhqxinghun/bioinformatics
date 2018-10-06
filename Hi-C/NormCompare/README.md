１. Source code for six H-C normalization methods. The code has been tested under Ubuntu14.04 LTS.

1> normSCNfun.R<br>
 Function:<br>
 　　　　Implementation of SCN method<br> 
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　./utilities/ scriptforSCN.sh (provided)<br>
 　　　　./SCN_sumV2.m (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　R.matlab (R package, )<br>
 　　　　MATALB (>= 8.6.0.267246 (R2015b))<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normSCNfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirtemp <- ./data/temp
         normSCNfun(dirin=dirin1, dirout=dirout1, dirtemp=dirtemp, species="hg19", chr="chr1", resolution=1000000)

2> normHiCNormfun.R<br>
 Function:<br>
 　　　　Implementation of HiCNorm method<br> 
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　BSgenome.Hsapiens.UCSC.hg19 (R package, 1.4.0)<br>
 　　　　BSgenome.Mmusculus.UCSC.mm9 (R package, 1.4.0)<br>
 　　　　HiTC (R package, 1.24.0)<br>
 　　　　Mappability files corresponding to hg19 and mm9 (Placed in extdata folder of HiTC)<br>
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

3> normICEfun.R<br>
 Function:<br>
 　　　　Implementation of ICE method<br>
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　HiTC (R package, 1.24.0)<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normICEfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirtemp <- ./data/temp
         normICEfun(dirin=dirin1, dirout=dirout1, dirtemp=dirtemp, species="hg19", chr="chr14", resolution=1000000)

4> normKRfun.R<br>
 Function:<br>
 　　　　Implementation of KR method<br>
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　./zerofilter.R (provided)<br>
 　　　　./zerorecover.R (provided)<br>
 　　　　./KR_bnewt.m (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　R.matlab (R package, )<br>
 　　　　MATALB (>= 8.6.0.267246 (R2015b))<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("zerofilter.R")
         source("zerorecover.R")
         source("normKRfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirtemp <- ./data/temp
         normKRfun(dirin=dirin1, dirout=dirout1, dirtemp=dirtemp, species="hg19", chr="chr1", resolution=1000000)

5> normchromoRfun.R<br>
 Function:<br>
 　　　　Implementation of chromoR method<br>
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　chromoR (R package, 1.0)<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normchromoRfun.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         normchromoRfun(dirin=dirin1, dirout=dirout1, species="hg19", chr="chr14", resolution=1000000)

6> normmultiHiCcomparefun.R<br>
 Function:<br>
 　　　　Implementation of multiHiCcompare method<br>
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　./gethicobjlist.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　multiHiCcompare (R package, 0.99.9)<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("gethicobjlist.R")
         source("normultiHiCcompare.R")
         dirin1 <- ./data/input/GM12878/chr1/1M
         dirin2 <- ./data/input/IMR90/chr1/1M
         dirin3 <- ./data/input/K562/chr1/1M
         dirout1 <- ./data/output/GM12878/chr1/1M
         dirout2 <- ./data/output/IMR90/chr1/1M
         dirout3 <- ./data/output/K562/chr1/1M
         normultiHiCcomparefun(dirin1=dirin1, dirin2=dirin2, dirin3=dirin3, dirout1=dirout1, dirout2=dirout2, dirout3=dirout3, species="hg19", chr="chr1", resolution=1000000)


２. Source code for logCPM transformation. The code has been tested under Ubuntu14.04 LTS.

 1> logCPM.R<br>
 Function:<br>
 　　　　Implementation of logCPM transformation<br>
 Depends:<br>
 　　　　./sparse2matrix.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 Example:
 
         rm(list=ls())
         source("sparse2matrix.R")
         source("logCPM.R")
         dirin1 <- ./data/output/GM12878/chr1/1M
         dirout1 <- ./data/logCPM/GM12878/chr1/1M
         logCPM(dirin=dirin1, dirout=dirout1, species="hg19", chr="chr1", resolution=1000000, bsparse=TRUE, bmodified=TRUE)


３. Datasets for examples. The data has been tested under Ubuntu14.04 LTS.<br>
 1> input<br>
 ./input/GM12878/chr1/1M<br>
 ./input/IMR90/chr1/1M<br>
 ./input/K562/chr1/1M<br>
2> output<br>
 ./output /GM12878/chr1/1M<br>
 ./output /IMR90/chr1/1M<br>
 ./output /K562/chr1/1M<br>
3> logCPM<br>
 ./logCPM /GM12878/chr1/1M<br>
 ./logCPM /IMR90/chr1/1M<br>
 ./logCPM /K562/chr1/1M<br>
