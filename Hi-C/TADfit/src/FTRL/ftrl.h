#ifndef FTRL_H
#define FTRL_H
#include "../datamanager.h"
#include "GSL/gsl_matrix_double.h"
#include <vector>

void ftrlalloc(gsl_matrix *&pftrlD, gsl_matrix *&pftrlY, gsl_matrix *&pftrlX, gsl_matrix *&pftrlB, gsl_matrix *&pftrlZ, gsl_matrix *&pftrlN, gsl_matrix *&pftrlP, int Ddim, int Ydim, int Xdim, int Bdim);
void ftrlfree(gsl_matrix *&pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP);
void ftrlgetsample(std::vector<contactmatrix> &sourdataset, int &curow, int &curcol, std::vector<contactdomain> &hiercd, gsl_matrix *pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, int &selcdnum);
void ftrlpredict(gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP, double &alpha, double &beta, double &l1, double &l2);
void ftrlgradient(gsl_matrix *pftrlB, int &curbrow, int &curbcol, double &curxvalue, double &curyvalue, double &curpvalue, double &lambda, double &gradient);
void ftrlupdate(gsl_matrix *pftrlY, gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP, std::vector<int> &contoidx, double &lambda, double &alpha);
void ftrloutput(std::vector<contactmatrix> &sourdataset, std::vector<contactdomain> &hiercd, std::vector<double> &regevar2, int selcdnum, gsl_matrix *pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, std::vector<int> &contoidx);

#endif
