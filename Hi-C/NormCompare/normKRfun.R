#function name:normKRfun
#function:read data from disks and use the method KR to normalize data
#parameters:
#             dirin:the path of input data
#             dirout:the path of output result
#             dirtemp:the path of temporary data
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
normKRfun <- function(dirin, dirout, dirtemp, species, chr, resolution) 
{
  require(R.matlab)
  #read data from disks and transform data to the data format of KR
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="KR")
  #function:use the method of KR to normalize data
  KRformat <- function(i)
  {
    filtedobj <- zerofilter(as.matrix(hicobjlist[[i]]))
    writeMat(paste0(dirtemp,"/","tempmatrix-preKR.mat"), A=filtedobj$matrix)
    system("./utilities/scriptforKR.sh -f ./KR_bnewt.m")
    hicnorm <- readMat(paste0(dirtemp,"/","tempmatrix-aftKR.mat"))
    recoveredhicnorm <- zerorecover(hicnorm$normatrix, filtedobj)
    write.table(round(recoveredhicnorm, digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-KR.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), KRformat)
}
