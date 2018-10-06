normultiHiCcomparefun <- function(dirin1, dirin2, dirin3, dirout1, dirout2, dirout3, species, chr, resolution)
{
  require(HiCcompare2)
  
  hicobjlist1 <- gethicobjlist(dirin1, species, chr, resolution, normethod="HiCcompare2")
  hicobjlist2 <- gethicobjlist(dirin2, species, chr, resolution, normethod="HiCcompare2")
  hicobjlist3 <- gethicobjlist(dirin3, species, chr, resolution, normethod="HiCcompare2")
  dim <- hicobjlist1$dim
  groups=c(rep(1,length(hicobjlist1$tablelist)),rep(2,length(hicobjlist2$tablelist)),rep(3,length(hicobjlist3$tablelist)))
  hicobjlist <- c(hicobjlist1$tablelist,hicobjlist2$tablelist,hicobjlist3$tablelist)
  rm(hicobjlist1, hicobjlist2, hicobjlist3)
  hicexp <- make_hicexp(data_list=hicobjlist, groups=groups, filter=FALSE)
  # hicexp <- cyclic_loess(hicexp, Plot = FALSE)
  hicexp <- fastlo(hicexp, Plot = FALSE)
  for(i in 1:length(hicobjlist))
  {
    if(groups[i] == 1)
      dirout <- dirout1
    else if(groups[i] == 2)
      dirout <- dirout2
    else if(groups[i] == 3)
      dirout <- dirout3
    hicsparsenorm <- data.frame(subset(hicexp@hic_table,select=c(2,3,i+4)))
    hicnorm <- sparse2matrix(hicsparsenorm,dim,resolution)
    write.table(round(hicnorm, digits=6), paste0(dirout,"/",names(hicobjlist)[i],"-HiCcompare2.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}