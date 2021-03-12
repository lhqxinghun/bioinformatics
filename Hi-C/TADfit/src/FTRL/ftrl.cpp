#include "ftrl.h"
#include "errlog.h"
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

void ftrlalloc(gsl_matrix *&pftrlD, gsl_matrix *&pftrlY, gsl_matrix *&pftrlX, gsl_matrix *&pftrlB, gsl_matrix *&pftrlZ, gsl_matrix *&pftrlN, gsl_matrix *&pftrlP, int Ddim, int Ydim, int Xdim, int Bdim)
{
    int i;
    pftrlD = gsl_matrix_alloc(1, Ddim);
    pftrlY = gsl_matrix_alloc(1, Ydim);
    pftrlX = gsl_matrix_alloc(1, Xdim);
    pftrlB = gsl_matrix_alloc(Xdim, Bdim);
    pftrlZ = gsl_matrix_alloc(Bdim, Xdim);
    pftrlN = gsl_matrix_alloc(Bdim, Xdim);
    pftrlP = gsl_matrix_alloc(1, Bdim);
    for(i = 0; i < Ddim; i++)
    {
        gsl_matrix_set(pftrlD, 0, i, log10(i+1)); //log(d+1)
        //gsl_matrix_set(pftrlD, 0, i, log10(i+0.3)); //log(d+0.5), log(d+0.3) 
    }
}

void ftrlfree(gsl_matrix *&pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP)
{
    gsl_matrix_free(pftrlD);
    gsl_matrix_free(pftrlY);
    gsl_matrix_free(pftrlX);
    gsl_matrix_free(pftrlB);
    gsl_matrix_free(pftrlZ);
    gsl_matrix_free(pftrlN);
    gsl_matrix_free(pftrlP);
}

void ftrlgetsample(std::vector<contactmatrix> &sourdataset, int &curow, int &curcol, std::vector<contactdomain> &hiercd, gsl_matrix *pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, int &selcdnum)
{
    int i, width;
    double curvalue;

    width = sourdataset.front().width;
    for(i = 0; i < (int)sourdataset.size(); i++)//给定坐标了，但是遍历样本
    {
        curvalue = *(sourdataset.at(i).pvalue + curow * width + curcol);
        gsl_matrix_set(pftrlY, 0, i, curvalue);//Y赋值
    }
    for(i = 0; i < selcdnum; i++)//X赋值需要遍历所有候选TADs，判断是否在范围内。
    {
        if((curow >= hiercd.at(i).top && curow <= hiercd.at(i).bottom && curcol >= hiercd.at(i).left && curcol <= hiercd.at(i).right) ||
           (curcol >= hiercd.at(i).top && curcol <= hiercd.at(i).bottom && curow >= hiercd.at(i).left && curow <= hiercd.at(i).right))
           //gsl_matrix_set(pftrlX, 0, i, log10(width-0.7)-gsl_matrix_get(pftrlD, 0, fabs(curow-curcol))); //线性衰减log(d+0.5)
            gsl_matrix_set(pftrlX, 0, i, log10(width)-gsl_matrix_get(pftrlD, 0, fabs(curow-curcol))); //线性衰减log(d+1)
            //gsl_matrix_set(pftrlX, 0, i, gsl_matrix_get(pftrlD, 0, width-fabs(curow-curcol)-1)); //log衰减
	    //gsl_matrix_set(pftrlX, 0, i, 1); //不衰减
        else
            gsl_matrix_set(pftrlX, 0, i, 0);
    }
}

void ftrlpredict(gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP, double &alpha, double &beta, double &l1, double &l2)
{
    int i, j, Xdim, Bdim, sign;
    double curxvalue, curzvalue, curnvalue, curbvalue;

    Bdim = pftrlB->size2; //context个数
    Xdim = pftrlX->size2; //候选TAD个数
    gsl_matrix_set_zero(pftrlP); //P矩阵置0
    for(i = 0; i < Xdim; i++) //遍历候选TAD
    {
        curxvalue = gsl_matrix_get(pftrlX, 0, i);
        if(curxvalue != 0)
        {
            for(j = 0; j < Bdim; j++)
            {
                curzvalue = gsl_matrix_get(pftrlZ, j, i);
                sign = curzvalue > 0? 1: -1;
                if((l1 > 0) && ((sign * curzvalue) <= l1))
                    gsl_matrix_set(pftrlB, i, j, 0);
                else
                {
                    curnvalue = gsl_matrix_get(pftrlN, j, i);
                    curbvalue = (sign * l1 - curzvalue) / ((beta + sqrt(curnvalue) / alpha + l2));
                    gsl_matrix_set(pftrlB, i, j, curbvalue);
                }
                gsl_matrix_set(pftrlP, 0, j, gsl_matrix_get(pftrlP, 0, j) + fabs(gsl_matrix_get(pftrlB, i, j)) * curxvalue);
            }
        }
    }
}

void ftrlgradient(gsl_matrix *pftrlB, int &curbrow, int &curbcol, double &curxvalue, double &curyvalue, double &curpvalue, double &lambda, double &gradient)
{
    int Bdim;
    double gradient1, gradient2, gradient3, curbvalue, leftbvalue, rightbvalue;;

    Bdim = pftrlB->size2;
    gradient1 = (curpvalue - curyvalue) * curxvalue;
    gradient2 = 0;
    gradient3 = 0;
    curbvalue = gsl_matrix_get(pftrlB, curbrow, curbcol);
    if(curbcol - 1 >= 0 && curbcol - 1 < Bdim)
    {
        leftbvalue = gsl_matrix_get(pftrlB, curbrow, curbcol - 1);
        gradient2 = fabs(curbvalue) > fabs(leftbvalue)? 1: -1;
    }
    if(curbcol + 1 >= 0 && curbcol + 1 < Bdim)
    {
        rightbvalue = gsl_matrix_get(pftrlB, curbrow, curbcol + 1);
        gradient3 = fabs(curbvalue) > fabs(rightbvalue)? 1: -1;
    }
    gradient = gradient1 + lambda * (gradient2 + gradient3);
    gradient = curbvalue > 0? gradient: (-1 * gradient);
}

void ftrlupdate(gsl_matrix *pftrlY, gsl_matrix *pftrlX, gsl_matrix *pftrlB, gsl_matrix *pftrlZ, gsl_matrix *pftrlN, gsl_matrix *pftrlP, std::vector<int> &contoidx, double &lambda, double &alpha)
{
    int i, j, Xdim, Ydim, conidx;
    double curxvalue, curyvalue, curpvalue, gradient, delta, curnvalue;

    Ydim = pftrlY->size2;
    Xdim = pftrlX->size2;
    for(i = 0; i < Xdim; i++)
    {
        curxvalue = gsl_matrix_get(pftrlX, 0, i);
        if(curxvalue != 0)
        {
            for(j = 0; j < Ydim; j++)
            {
                conidx = contoidx.at(j);
                curyvalue = gsl_matrix_get(pftrlY, 0, j);
                curpvalue = gsl_matrix_get(pftrlP, 0, conidx);
                curnvalue = gsl_matrix_get(pftrlN, conidx, i);
                ftrlgradient(pftrlB, i, conidx, curxvalue, curyvalue, curpvalue, lambda, gradient);
                delta = (sqrt(curnvalue + pow(gradient, 2)) - sqrt(curnvalue)) / alpha;
                gsl_matrix_set(pftrlZ, conidx, i, gsl_matrix_get(pftrlZ, conidx, i) + gradient - delta * gsl_matrix_get(pftrlB, i, conidx));
                gsl_matrix_set(pftrlN, conidx, i, curnvalue + pow(gradient, 2));
            }
        }
    }
}

void ftrloutput(std::vector<contactmatrix> &sourdataset, std::vector<contactdomain> &hiercd, std::vector<double> &regevar2, int selcdnum, gsl_matrix *pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, std::vector<int> &contoidx)
{
    int i, j, m, n, k, width, height, Ydim, Xdim, Bdim, repnum;
    double sumyvalue, sumbvalue, curpvalue, curxvalue, curyvalue, rsquared;
    double *paveyvalue, *pavebvalue, *presvalue, *ptotvalue;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    Ydim = pftrlY->size2;//样本个数
    Xdim = pftrlX->size2;//候选TAD个数
    Bdim = hiercd.front().conum;//context个数
    repnum = hiercd.front().repnum;//重复次数
    paveyvalue = new (std::nothrow) double[Ydim];
    if(paveyvalue == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    pavebvalue = new (std::nothrow) double[Xdim * Bdim];
    if(pavebvalue == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    presvalue = new (std::nothrow) double[Ydim];
    if(presvalue == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    memset(presvalue, 0, sizeof(double) * Ydim);
    ptotvalue = new (std::nothrow) double[Ydim];
    if(ptotvalue == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    memset(ptotvalue, 0, sizeof(double) * Ydim);
    //Calculate R squared for each dataset
    for(i = 0; i < Ydim; i++)//遍历每个样本，结果是每个样本都有了一个平均值，保存地址paveyvalue
    {
        sumyvalue = 0;
        for(j = 0; j < height; j++)
        {
            for(m = 0; m<= j; m++)//行，列遍历
                sumyvalue += *(sourdataset.at(i).pvalue + j * width + m);
        }
        *(paveyvalue + i) = sumyvalue / ((width + 1) * height / 2);
    }
    for(i = 0; i < Xdim; i++)//遍历每个候选TAD
    {
        for(j = 0; j < Bdim; j++)//遍历每个context，结果是每个候选TAD的每个context都具有了一个平均beta值，保存地址pavebvalue
        {
            sumbvalue = 0;
            for(m = 0; m < repnum; m++)//重复次数遍历
            {
                sumbvalue += *(hiercd.at(i).pvalue + j * repnum + m);
            }
            *(pavebvalue + j * Xdim + i) = sumbvalue / repnum;
        }
    }
    for(i = 0; i < Bdim; i++)//遍历每个contex
    {
        for(j = 0; j < height; j++)
        {
            for(m = 0; m <= j; m++)//行列遍历
            {
                curpvalue = 0;
                ftrlgetsample(sourdataset, j, m, hiercd, pftrlD, pftrlY, pftrlX, selcdnum);
                for(n = 0; n < selcdnum; n++)
                {
                    curxvalue = gsl_matrix_get(pftrlX, 0, n);
                    if(curxvalue != 0)
                    {
                        curpvalue += *(pavebvalue + i * Xdim + n) * curxvalue;
                    }
                }//计算得到当前点的pvalue(预测值）,根据当前context，当前IF点的X向量和beta的平均值。
                for(k = 0; k < Ydim; k++)//遍历所有样本
                {
                    if(contoidx.at(k) == i)//如果这个样本属于当前context
                    {
                        curyvalue = *(sourdataset.at(k).pvalue + j * width + m);//把这个点的Y值取出来
                        *(presvalue + k) += pow(curyvalue - curpvalue, 2);      //累加Y和P的残差平方和，存在presvalue里
                        *(ptotvalue + k) += pow(curyvalue - *(paveyvalue + k), 2);//累加得到总离差平方和，存在ptotvalue里
                        *(sourdataset.at(k).pvalue + j * width + m) = curpvalue;//完了以后会更新掉当前样本的IF值?
                        if(j != m)
                            *(sourdataset.at(k).pvalue + m * width + j) = curpvalue;//转置IF矩阵
                    }
                }
            }
        }
    }
    for(i = 0; i < Ydim; i++)
    {
        rsquared = 1- (*(presvalue + i) / *(ptotvalue + i));
        regevar2.push_back(rsquared);
    }
    delete []paveyvalue;
    delete []pavebvalue;
    delete []presvalue;
    delete []ptotvalue;
}