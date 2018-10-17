#function name:normultiHiCcomparefun
#function:read data from disks and use the method multiHiCcompare to normalize data
#parameters:
#             dirin1:the path of the first input data
#             dirin2:the path of the second input data
#             dirin3:the path of the third input data  
#             dirout1:the path of the first output result
#             dirout2:the path of the second output result
#             dirout3:the path of the third output result
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
normultiHiCcomparefun <- function(dirin1, dirin2, dirin3, dirout1, dirout2, dirout3, species, chr, resolution)
{
  require(HiCcompare2)
  #read data from disks and transform data to the data format of multiHiCcompare
  hicobjlist1 <- gethicobjlist(dirin1, species, chr, resolution, normethod="multiHiCcompare")
  hicobjlist2 <- gethicobjlist(dirin2, species, chr, resolution, normethod="multiHiCcompare")
  hicobjlist3 <- gethicobjlist(dirin3, species, chr, resolution, normethod="multiHiCcompare")
  dim <- hicobjlist1$dim
  groups=c(rep(1,length(hicobjlist1$tablelist)),rep(2,length(hicobjlist2$tablelist)),rep(3,length(hicobjlist3$tablelist)))
  hicobjlist <- c(hicobjlist1$tablelist,hicobjlist2$tablelist,hicobjlist3$tablelist)
  rm(hicobjlist1, hicobjlist2, hicobjlist3)
  #function:use the method of multiHiCcompare to normalize data
  hicexp <- make_hicexp(data_list=hicobjlist, groups=groups, filter=FALSE)
  
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
    write.table(round(hicnorm, digits=6), paste0(dirout,"/",names(hicobjlist)[i],"-multiHiCcompare.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}
