#ifndef TOPDOM_H
#define TOPDOM_H

#include <vector>

//Cd type
typedef enum
{
  gap,
  domain,
  boundary
} cdtype;
//Region
typedef struct
{
  cdtype type;
  int start;
  int end;
  double stpvalue;
  double edpvalue;
} cdcps;

void getdiamondmean(double *pmatrix, int width, int height, int size, double *pdiamean);
void getgapindex(double *pmatrix, int width, int height, std::vector<int> &gapidx);
void getproregion(std::vector<int> &gapidx, int width, int height, int minsize, std::vector<cdcps> &proregion);
void diameannorm(double *psourdata, int start, int end, double *pdestdata);
void getchangepoint(double *pnordiamean, int start, int end, std::vector<int> &changepoint);
void detectextreme(double *pdiamean, int len, std::vector<cdcps> &proregion, double *plocalext);
double WilcoxRanksumTest(double *pvector1, int size1, double *pvector2, int size2);
void getpvalue(double *pmatrix, int width, int height, std::vector<cdcps> &proregion, int size, double *pvalue);
void bintodomain(double *plocalext, int width, int height, double *pvalue, bool bstatFilter, std::vector<cdcps> &diagcd);
void getdiagchangepoint(double *pmatrix, int width, int height, int winsize, bool bstatFilter, std::vector<cdcps> &diagcd);

#endif
