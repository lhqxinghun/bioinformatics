# TADfit

TADfit is a multivariate linear regression model for profiling of hierarchical chromatin domains, which tries to fit the interaction frequencies in Hi-C contact matrix with and without replicates using all-possible hierarchical TADs, and the regression coefficient for each TAD is obtained by an online learning solver called Follow-The-Regularized -Leader (FTRL) to determine the significant ones from them.

## Installation
 The TADfit was built on R 3.6.2 and compiled with gcc 7.4.0, and the same version of platforms is recommended for users.
 1) If necessary, install the required dependencies in R. 
```
install.packages(c('Matrix', 
                'perm', 
                'matrixStats', 
                'Rcpp', 
                'HiCcompare', 
                'plyr', 
                'zoo'))
```
 2) Install TADfit from source package

  (Note that gcc with corresponding version are needed)
```
  install.packages("path to TADfit_1.1.tar.gz", repos = NULL, type = "source") 
```  
  
 3) Or install TADfit from compiled package
```
#linux
install.packages('path to TADfit_1.1_R_x86_64-pc-linux-gnu.tar.gz')
```
```
#Windows
install.packages('path to TADfit_1.1.zip', repos = NULL, type = "win.binary")
```

## Usage
The hierarchical TADs on replicate Hi-C data are called mainly by the  function ==calltad== in TADfit. Other auxiliary functions can be accessable in the package too.
### Input data
Before runing the TADfit, Users need to prepare the input data. The input data to TADfit is a group of replicate normalized Hi-C contact matrices under the same scale. Here we recommend iterative correction and eigenvector decomposition (ICE) followed by Counts per million (CPM) for normalization. These matrix files need to be put on the same directory.

### Parameters
- **dirin**: path to a directory which contains replicate Hi-Ccontact matrix file(s)
- **resolution**: bin size of the H-C contact matrix
- **start**: start coordinate of the regionof chromosome
- **end**: start coordinate of the regionof chromosome
- **specied**: "hg19" for human or "mm9"for mouse
- **chr**: chromosome number, e.g. "chr1"
- **winsize**:window size when calling TAD boundaries. This parameter is from TopDom (default 5)
- **minsizescale**: a scale factor indicates the minimum size of a hierarchical TAD (default 0.01)
- **maxsizescale**: a scale factor indicates the maximum size of a hierarchical TAD (default 0.5)
- **repeatnum**: the number of repeat when solving the mode l(default 10)
- **iteration**: the number of iteration when solving the model (default 5)
- **alpha, beta**: parameters for pre-coordinating learning rate when optimizing the model (default 0.01 for alpha and 1 for beta)
- **l**: adjust the strength of l1 regularization (defalut 2)

### Output
The model(calltad function) outputs a dataframe with four colomns
- **start**: start coordinate (in bin) of a hierarchical TAD
- **end**: end coordinate (in bin)of a hierarchical TAD
- **b_avg**: the regression coefficients estimated by FTRL
- **pvalue**: the p-value for permutation test

An example of output:
```
start   end     b_avg       pvalue
0       7       0.11833     0.00082
2       7       0.14427     0.00079
7       18      0.22934     0.00054
7       33      0.03613     0.00157
...
```

### Example
Call hierarchical TADs from the sample Hi-C data
```
library(TADfit)
dirin <- "/home/biology/data/expdata/GM12878-chr1"
resolution <- 25000
start <- 0
end <- 5000000
TADs <- calltad(dirin, resolution, start, end)  #Other parameters are default
```
## Reference
xxx


## Contact
xxx

