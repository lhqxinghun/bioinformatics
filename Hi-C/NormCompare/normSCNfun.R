normSCNfun <- function(dirin, dirout, dirtemp, species, chr, resolution) 
{
  require(R.matlab)
  
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="SCN")
  
  SCNformat <- function(i)
  {
    writeMat(paste0(dirtemp,"/","tempmatrix-preSCN.mat"), DataF=hicobjlist[[i]])
    system("./utilities/scriptforSCN.sh -f ./SCN_sumV2.m")
    hicnorm <- readMat(paste0(dirtemp,"/","tempmatrix-aftSCN.mat"))
    write.table(round(as.matrix(hicnorm$DataFN), digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-SCN.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), SCNformat)
}