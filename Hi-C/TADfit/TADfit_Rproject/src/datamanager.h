#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include "FTRL/GSL/gsl_matrix_double.h"
#include <string>
#include <vector>

//Lenght of string for context
const int ctlen = 101;
//Size for memeory allocation
const int g_cmemory = 20000;
const int g_cmemstep = 10000;

//Geometrical Type of contact matrix
typedef enum
{
  rectangle,
  uppertri,
  lowertri,
  diagonal
} geotype;

//Methods for pretreatment of contact matrix
typedef enum
{
    log2m,
    log10m,
    medfm,
    znorm,
    dnorm
} pretrementmethod;

//Contact matrix
typedef struct
{
  char context[ctlen];
  geotype type;
  int width;
  int height;
  double *pvalue;
} contactmatrix;

//potential change point
typedef struct
{
  int diagpos;
  double score;
} changepoint;

//Potential contact domain
typedef struct
{
  int left;
  int right;
  int top;
  int bottom;
  int conum;
  int repnum;
  double *pvalue;
} contactdomain;

//Parameter for median filtering and distance normalization
typedef struct
{
    bool bmedifilter;
    bool bicnorm;
    int maxiter;
} paranorm;

//Parameters for TopDom and TADs filter
typedef struct
{
    int winsize;
    bool bstatfilter;
    bool bTADsfilter;
    double minsizescale;
    double maxsizescale;
} parahierTADs;

//Parameters for FTRL regression
typedef struct
{
    int repeatnum;
    int iteration;
    double lambda;
    double alpha;
    double beta;
    double l1;
    double l2;
} paraFTRLreg;

class DataManager
{
public:
    DataManager();
    ~DataManager();

public:
    void cmpretrement(std::vector<contactmatrix> &sourdataset, pretrementmethod selmethod);
    void icnormalization(std::vector<contactmatrix> &sourdataset, int maxiter);
    void findiagcpbyTopDom(std::vector<contactmatrix> &sourdataset, int winsize, bool bstatFilter, std::vector<changepoint> &diagcp);
    void diagcptohiercd(std::vector<contactmatrix> &sourdataset, std::vector<changepoint> &diagcp, int repeatnum, bool bTADsfilter, double minsizescale, double maxsizescale, std::vector<contactdomain> &hiercd);
    void hiercdfilter(std::vector<contactmatrix> &sourdataset, double minsizescale, double maxsizescale, std::vector<contactdomain> &hiercd);
    void ftrlsolver(std::vector<contactmatrix> &sourdataset, std::vector<contactdomain> &hiercd, std::vector<double> &regevar2,
                    int repeatnum, int iteration, double lambda, double alpha, double beta, double l1, double l2);

    //temporary function
    inline void linetovector(std::string &line, std::vector<std::string> &fields, char delimiter);
    int readtable(double *&pdestdata, int &width, int &height, std::string sourfilepath, bool bnewmem);
    void writetable(std::string destfilepath, double *psourdata, int width, int height);

private:
    void cmlog2(double *psourdata, int width, int height, double *pdestdata);
    void cmlog10(double *psourdata, int width, int height, double *pdestdata);
    void cmedianfilter(double *psourdata, int width, int height, double *pdestdata);
    void cmznormalization(double *psourdata, int width, int height, double *pdestdata);
    void sdnormalization(std::vector<contactmatrix> &sourdataset);
    void calpseudoreference(std::vector<contactmatrix> &sourdataset, contactmatrix &pserefdataset);
    void mapcontextindex(std::vector<contactmatrix> &sourdataset, std::vector<int> &contextoindex, std::vector<std::string> &indextocontext);
};

#endif // DATAMANAGER_H
