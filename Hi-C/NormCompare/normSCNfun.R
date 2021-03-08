#function name:normSCNfun
#function:read data from disks and use the method SCN to normalize data
#parameters:
#             dirin:the path of input data
#             dirout:the path of output result
#             dirtemp:the path of temporary data
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
normSCNfun <- function(dirin, dirout, dirtemp, species, chr, resolution) 
{
  require(R.matlab)
  #read data from disks and transform data to the data format of SCN
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="SCN")
  #function:use the method of SCN to normalize data
  SCNformat <- function(i)
  {
    writeMat(paste0(dirtemp,"/","tempmatrix-preSCN.mat"), DataF=hicobjlist[[i]])
    system("./utilities/scriptforSCN.sh -f ./SCN_sumV2.m")
    hicnorm <- readMat(paste0(dirtemp,"/","tempmatrix-aftSCN.mat"))
    write.table(round(as.matrix(hicnorm$DataFN), digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-SCN.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), SCNformat)
}
