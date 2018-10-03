normKRfun <- function(dirin, dirout, dirtemp, species, chr, resolution) 
{
  require(R.matlab)
  
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="KR")
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