lower.tri <- function(filein.y, filein.p, repnum, selycon, selyrep, selpcon, conum, selpiter, fileout)
{
  library(gplots)
  require(pheatmap)
  require(data.table)  

  yvalue <- as.matrix(fread(filein.y, header=F))
  pvalue <- as.matrix(fread(filein.p, header=F))
  selyrow <- (selycon-1)*repnum+selyrep
  selprow <- (selpiter-1)*conum+selpcon
  #print(dim(yvalue))
  #print(dim(pvalue))
  totalnum<-dim(yvalue)[2]
  #print(totalnum)
  originalval <- as.matrix(yvalue[selyrow,])
  fittedval <- as.matrix(pvalue[selprow,])
  total<- nrow(fittedval)
  f<- function(x,a,b,c) a*x^2+b*x+c
  a<- 1
  b<- 1
  c<- (-2)*totalnum
  num<- uniroot(f, c(0,10000),a=a, b=b, c=c, tol=0.00001)
  num<- num$root
  print(num)
  tri<- matrix(0,nrow=num,ncol=num)
  n<- 1
  for (i in 1:num)
  {
    for (j in 1:i)
    {
      tri[i,j]<-fittedval[n,1]
      n<- n+1
    }
  }
  n<-1
  #tri<- tri+t(tri)
  for(j in 1:num){
    for(i in 1:j){
      tri[i,j]<-originalval[n,1]
      n<- n+1
    }
  }
  for (i in 1:num) 
  {
  tri[i,i]<-tri[i,i]/2
  }

  #tri[tri>1.8]=1.8
  #tri <- log10(tri+1)
  fileout.heatmap <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,"-iter",selpiter, ".TIFF")
  tiff(fileout.heatmap, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=400)
  par(mar=c(9, 11, 0.2, 0), mgp=c(6,3,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
  breaks = seq(min(tri),  max(tri), length.out=129)
  pheatmap(tri,breaks=breaks,col=colorpanel(128,"lightyellow","red"),border_color=NA,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,fontsize=4)
  #pheatmap(log10(tri+1),col=colorpanel(128,"lightyellow","red"),border_color=NA,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,fontsize=4)
  dev.off()
  
  # tiff(fileoutlist[[i]], width = 2.46, height = 2.46, units = "cm", pointsize=1, res=350)
  # #tiff(fileoutlist[[i]], width = 6, height = 6, units = "cm", pointsize=1, res=350)
  # coordstr <-  sprintf("%.2fM", seq(start/1000000, end/1000000, length.out=nrow(submatrixIFlist[[i]])))
  # breaks = seq(minlogIF,  maxlogIF, length.out=129)
  # pheatmap(submatrixIFlist[[i]],breaks=breaks,col=colorpanel(128,"lightyellow","red"),border_color=NA,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,fontsize=4)
  # #pheatmap(submatrixIFlist[[i]],breaks=breaks,col=colorpanel(128,"lightyellow","red"),border_color=NA,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,legend=FALSE)
  # dev.off()
}

