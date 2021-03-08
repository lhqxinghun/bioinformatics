#function name:gethicobjlist
#function:read data from disk and transform data to the certain format based on the given parameters,such as the resolution,the normethod and so on
#parameters:
#           dir:the path of input data
#           species:hg19 or mm9
#           chr:the numberID of the chromosome 
#           resolution:the resolution of input data
#           normethod:the method of normalization, including SCN, HiCNorm, ICE, KR, chromoR and multiHiCcompare
gethicobjlist <- function(dir, species, chr, resolution, normethod) 
{
  #hg19 or mm9
  switch(species,
    hg19=
    {
      chrlens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                  146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                  102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                  155270560, 59373566)
      names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                          "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                          "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                          "chrX", "chrY")
      dim <- floor(chrlens[chr] / resolution) + 1      
    },
    mm9=
    {
      chrlens = c(197195432, 181748087, 159599783, 155630120, 152537259, 149517037,	152524553,
                  131738871,	124076172,	129993255,	121843856, 121257530,	120284312,	125194864,
                  103494974,	98319150,	95272651,	90772031,	61342430,	166650296,	15902555,	16300)
      names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                          "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                          "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM")
      dim <- floor(chrlens[chr] / resolution) + 1      
    }    
  )
  #read data from disks and transform data to matrixs
  tablelist <- sapply(list.files(dir, pattern=".txt", full.names=TRUE), read.table, simplify = FALSE, USE.NAMES = TRUE)
  matrixlist <- lapply(tablelist, sparse2matrix, dim, resolution)
  namesplit <- strsplit(basename(names(tablelist)),split=".",fixed=TRUE)
  names(tablelist) <- unlist(lapply(namesplit,head,1))
  namesplit <- strsplit(basename(names(matrixlist)),split=".",fixed=TRUE)
  names(matrixlist) <- unlist(lapply(namesplit,head,1))
  
  #transform data to the certain format based on the given normethod
  switch(normethod,
    SCN=
    {
      options(scipen = 100)
      hicobjlist <- matrixlist
    },
    HiCNorm=
    {
      options(scipen = 100)
      start <- seq(0,chrlens[chr],resolution)
      end <- c(start[1:(length(start)-1)]+resolution-1, chrlens[chr])
      names <- paste0("HIC_bin", 1:dim, "|hg19mm9|", chr, ":", start, "-", end)
      setrowcolnames <- function(x)
      {
        rownames(x) <- names
        colnames(x) <- names
        return(x)
      }
      hicobjlist <- lapply(matrixlist, setrowcolnames)
    },
    ICE=
    {
      options(scipen = 100)
      start <- seq(0,chrlens[chr],resolution)
      end <- c(start[1:(length(start)-1)]+resolution-1, chrlens[chr])
      names <- paste0("HIC_bin", 1:dim, "|hg19mm9|", chr, ":", start, "-", end)
      setrowcolnames <- function(x)
      {
        rownames(x) <- names
        colnames(x) <- names
        return(x)
      }
      hicobjlist <- lapply(matrixlist, setrowcolnames)
    },
    KR=
    {
      options(scipen = 100)
      hicobjlist <- matrixlist
    },
    chromoR=
    {
      options(scipen = 100)
      start <- seq(0,chrlens[chr],resolution)
      end <- c(start[1:(length(start)-1)]+resolution-1, chrlens[chr])
      seg <- data.frame(chr=chr, start=start, end=end)
      hicobjlist <- list(matrixlist=matrixlist, seg=seg)
    },
    multiHiCcompare=
    {
      options(scipen = 100)
      tablelist <- lapply(tablelist, function(x) cbind(substr(chr,4,nchar(chr)), x))
      hicobjlist <- list(tablelist=tablelist, dim=dim)
    }
  )
  return (hicobjlist)
}
