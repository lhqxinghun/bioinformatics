#include "TopDom.h"
#include "errlog.h"
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <float.h>

void getdiamondmean(double *pmatrix, int width, int height, int size, double *pdiamean)
{
    int i, m, n, lowerbound, upperbound;
    double sumvalue;

    for(i = 0; i < height; i++)
    {
        if(i != height -1)
        {
            lowerbound = (i - size + 1) > 0? (i - size + 1) : 0;
            upperbound = (i + size) < (height - 1)? (i + size) : (height - 1);
            sumvalue = 0;
            for(m = lowerbound; m <= i; m++)
            {
                for(n = i + 1; n <= upperbound; n++)
                    sumvalue += *(pmatrix + m * width + n);
            }
            *(pdiamean + i) = sumvalue / ((i - lowerbound + 1) * (upperbound - i));
        }
        else
            *(pdiamean + i) = 0;
    }
}

void getgapindex(double *pmatrix, int width, int height, std::vector<int> &gapidx)
{
    int i, j, m , n;
    double sumvalue;

    i = 0;
    while(i < height - 1)
    {
        for(j = i + 1; j < height; j++)
        {
            sumvalue = 0;
            for(m = i; m <= j; m++)
                for(n = i; n <= j; n++)
                    sumvalue +=  *(pmatrix + m * width + n);
            if(sumvalue == 0)
            {
                for(m = i; m <= j; m++)
                    gapidx.push_back(m);
            }
            else
                break;
        }
        i = j;
    }
}

void getproregion(std::vector<int> &gapidx, int width, int height, int minsize, std::vector<cdcps> &proregion)
{
    int i, j, start;
    cdcps tmpregion;
    std::vector<int> proset;
    std::vector<cdcps> tmproregion;

    if(width == height)
    {
        for(i = 0; i < height; i++)
        {
            if(std::find(gapidx.begin(), gapidx.end(), i) == gapidx.end())
               proset.push_back(i);
        }
        i = 0;
        while(i < (int)proset.size() - 1)
        {
            start = proset[i];
            for(j = i + 1; j < (int)proset.size(); j++)
            {
                if(proset[j] - proset[j-1] <= 1)
                    continue;
                else
                {
                    tmpregion.start = start;
                    tmpregion.end =  proset[j-1];
                    tmproregion.push_back(tmpregion);
                    i = j;
                    break;
                }
            }
            if(j >= (int)proset.size() - 1)
            {
                tmpregion.start = start;
                tmpregion.end =  proset[j-1];
                tmproregion.push_back(tmpregion);
                break;
            }
        }
        for(i = 0; i < (int)tmproregion.size(); i++)
        {
            if(abs(tmproregion.at(i).end - tmproregion.at(i).start) >= minsize)
                    proregion.push_back(tmproregion.at(i));
        }
    }
}

void diameannorm(double *psourdata, int start, int end, double *pdestdata)
{
    int i, len;
    double sumvalue, scale;
    double *pcursourdata, *pcurdestdata, *pdiff;

    len = end - start + 1;
    pcursourdata = psourdata + start;
    pcurdestdata = pdestdata + start;
    sumvalue = 0;
    pdiff = new (std::nothrow) double[len - 1];
    if(pdiff == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    *pcurdestdata = *pcursourdata;
    for(i = 0; i < len - 1; i++)
    {
        *(pdiff + i) = *(pcursourdata + i + 1) - *(pcursourdata + i);
        sumvalue += fabs(*(pdiff + i));
    }
    scale = (len - 1) / sumvalue;
    for(i = 1; i < len; i++)
        *(pcurdestdata + i) = *(pcurdestdata + i - 1) + *(pdiff + i - 1) * scale;
    delete []pdiff;
}

void getchangepoint(double *pnordiamean, int start, int end, std::vector<int> &changepoint)
{
    int i, j, m, len;
    double sumvalue;
    double *pcurnordiamean, *pfv, *pev;

    len = end - start + 1;
    pcurnordiamean = pnordiamean + start;
    pfv = new (std::nothrow) double[len];
    if(pfv == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    else
    {
        for(i = 0; i < len; i++)
            *(pfv + i) = nan("");
    }
    pev = new (std::nothrow) double[len];
    if(pev == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    else
    {
        for(i = 0; i < len; i++)
            *(pev + i) = nan("");
    }
    changepoint.push_back(start);
    *pfv = 0;
    i = 0;
    while(i < len - 1)
    {
        j = i + 1;
        pfv[j] = sqrt(pow(j - i, 2) + pow(*(pcurnordiamean + j) - *(pcurnordiamean + i), 2));
        while(j < len - 1)
        {
            j++;
            sumvalue = 0;
            for(m = i + 1; m <= j - 1; m++)
                sumvalue += fabs((*(pcurnordiamean + j) - *(pcurnordiamean + i)) * m - (j - i) * *(pcurnordiamean + m) - (i * *(pcurnordiamean + j)) + (j * *(pcurnordiamean + i)));
            *(pev + j) = sumvalue / sqrt(pow(j - i, 2) + pow(*(pcurnordiamean + j) - *(pcurnordiamean + i), 2));
            *(pfv + j) = sqrt(pow(j - i, 2) + pow(*(pcurnordiamean + j) - *(pcurnordiamean + i), 2)) -  *(pev + j);
            if(std::isnan(*(pfv + j)) || std::isnan(*(pfv + j - 1)))
            {
                j--;
                changepoint.push_back(start + j);
                break;
            }
            if(*(pfv + j) < *(pfv + j - 1))
            {
              j--;
              changepoint.push_back(start + j);
              break;
            }
        }
        i = j;
    }
    changepoint.push_back(end);
    delete []pfv;
    delete []pev;
}

void detectextreme(double *pdiamean, int len, std::vector<cdcps> &proregion, double *plocalext)
{
    int i, j, m, start, end, minidx, maxidx;
    double minvalue, maxvalue, minval, maxval;
    double *pnordiamean;
    std::vector<int> changepoint;

    minidx = 0;
    maxidx = 0;
    pnordiamean = new (std::nothrow) double[len];
    if(pnordiamean == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    else
        memset(pnordiamean, 0, sizeof(double) * len);
    for(i = 0; i < (int)proregion.size(); i++)
    {
        start = proregion.at(i).start;
        end = proregion.at(i).end;
        changepoint.clear();
        memset(plocalext + start, 0, sizeof(double) * (end - start + 1));
        if(end - start + 1 <= 3)
        {
            minvalue = DBL_MAX;
            maxvalue = -DBL_MAX;
            for(j = start; j <= end; j++)
            {
                if(*(pdiamean + j) < minvalue)
                {
                    minvalue = *(pdiamean + j);
                    minidx = j;
                }
                if(*(pdiamean + j) > maxvalue)
                {
                    maxvalue = *(pdiamean + j);
                    maxidx = j;
                }
            }
            *(plocalext + minidx) = -1;
            *(plocalext + maxidx) = 1;
        }
        else
        {
            diameannorm(pdiamean, start, end, pnordiamean);
            getchangepoint(pnordiamean, start, end, changepoint);
            if((int)changepoint.size() > 2)
            {
                for(j = 1; j < (int)changepoint.size() - 1; j++)
                {
                    if(*(pnordiamean + changepoint.at(j)) >= *(pnordiamean + changepoint.at(j) - 1) && *(pnordiamean + changepoint.at(j)) >= *(pnordiamean + changepoint.at(j) + 1))
                        *(plocalext + changepoint.at(j)) = 1;
                    else if(*(pnordiamean + changepoint.at(j)) < *(pnordiamean + changepoint.at(j) - 1) && *(pnordiamean + changepoint.at(j)) < *(pnordiamean + changepoint.at(j) + 1))
                        *(plocalext + changepoint.at(j)) = -1;
                    minval = *(pnordiamean + changepoint.at(j-1)) < *(pnordiamean + changepoint.at(j))? *(pnordiamean + changepoint.at(j-1)): *(pnordiamean + changepoint.at(j));
                    maxval = *(pnordiamean + changepoint.at(j-1)) > *(pnordiamean + changepoint.at(j))? *(pnordiamean + changepoint.at(j-1)): *(pnordiamean + changepoint.at(j));
                    minvalue = DBL_MAX;
                    maxvalue = -DBL_MAX;
                    for(m = changepoint.at(j-1); m <= changepoint.at(j); m++)
                    {
                        if(*(pnordiamean + m) < minvalue)
                        {
                            minvalue = *(pnordiamean + m);
                            minidx = m;
                        }
                        if(*(pnordiamean + m) > maxvalue)
                        {
                            maxvalue = *(pnordiamean + m);
                            maxidx = m;
                        }
                    }
                    if(minvalue < minval)
                        *(plocalext + minidx) = -1;
                    if(maxvalue > maxval)
                        *(plocalext + maxidx) = 1;
                }
            }
        }
    }
    delete []pnordiamean;
}

//Returns the cumulative probability of x of a gaussian distribution
double normcdf(double x, double mu, double sigma)
{
    int sign;
    double t, y;

    //standardization
    x = (x - mu) / sigma;
    //Constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    //Save the sign of x
    sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    //A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    return 0.5*(1.0 + sign*y);
}

double WilcoxRanksumTest(double *pvector1, int size1, double *pvector2, int size2)
{
    bool tmpflag;
    bool *pflag;
    int i, j, m, n, tmpidx, totalnum, tieindex;
    int *ptiesize;
    double tmprank, R, U, mu, sigma, cdfvalue;
    double *prank;

    tmpidx = 0;
    for(i = 0; i < size1; i++)
    {
        if(std::isnan(*(pvector1 + i)))
            continue;
        else
        {
            *(pvector1 + tmpidx) = *(pvector1 + i);
            tmpidx++;
        }
    }
    size1 = tmpidx;
    tmpidx = 0;
    for(i = 0; i < size2; i++)
    {
        if(std::isnan(*(pvector2 + i)))
            continue;
        else
        {
            *(pvector2 + tmpidx) = *(pvector2 + i);
            tmpidx++;
        }
    }
    size2 = tmpidx;
    totalnum = size1 + size2;
    prank = new (std::nothrow) double[totalnum];
    if(prank == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    pflag = new (std::nothrow) bool[totalnum];
    if(pflag == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    ptiesize = new (std::nothrow) int[totalnum];
    if(pflag == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < size1; i++)
    {
        *(prank + i) = *(pvector1 + i);
        *(pflag + i) = false;
    }
    for(i = 0; i < size2; i++)
    {
        *(prank + i + size1) = *(pvector2 + i);
        *(pflag + i + size1) = true;
    }
    //sort
    if(totalnum != 1)
    {
        i = 2;
        do
        {
            n = i;
            while(n != 1)
            {
                m = n / 2;
                if(*(prank + m - 1) >= *(prank + n - 1))
                {
                    n = 1;
                }
                else
                {
                    tmprank = *(prank + m - 1);
                    *(prank + m - 1) = *(prank + n - 1);
                    *(prank + n - 1) = tmprank;
                    tmpflag = *(pflag + m - 1);
                    *(pflag + m - 1) = *(pflag + n - 1);
                    *(pflag + n - 1) = tmpflag;
                    n = m;
                }
            }
            i = i + 1;
        }
        while(i <= totalnum);
        i = totalnum - 1;
        do
        {
            tmprank = *(prank + i);
            *(prank + i) = *prank;
            *prank = tmprank;
            tmpflag =*(pflag + i);
            *(pflag + i) =*pflag;
            *pflag = tmpflag;
            n = 1;
            while(n != 0)
            {
                m = 2 * n;
                if(m > i)
                {
                    n = 0;
                }
                else
                {
                    if(m < i)
                    {
                        if(*(prank + m) > *(prank + m - 1))
                        {
                            m = m + 1;
                        }
                    }
                    if(*(prank + n - 1) >= *(prank + m - 1))
                    {
                        n = 0;
                    }
                    else
                    {
                        tmprank = *(prank + m - 1);
                        *(prank + m - 1) = *(prank + n - 1);
                        *(prank + n - 1) = tmprank;
                        tmpflag = *(pflag + m - 1);
                        *(pflag + m - 1) = *(pflag + n - 1);
                        *(pflag + n - 1) = tmpflag;
                        n = m;
                    }
                }
            }
            i = i - 1;
        }
        while(i >= 1);
    }
    //Compute tied ranks
    i = 0;
    tieindex = 0;
    while(i < totalnum)
    {
        j = i + 1;
        while(j <= totalnum - 1)
        {
            if(*(prank + j) != *(prank + i))
            {
                break;
            }
            j = j + 1;
        }
        for(m = i; m <= j - 1; m++)
        {
            *(prank + m) = 1 + (i + j - 1) / (double)2.0;
        }
        *(ptiesize + tieindex) = j - i;
        tieindex++;
        i = j;
    }
    //Compute R
    R = 0;
    for(i = 0; i < totalnum; i++)
    {
        if(*(pflag + i) == 0)
        {
            R += *(prank + i);
        }
    }
    //Compute U with continuity correction
    U = R - size1 * (size1 + 1) / (double)2.0 + 0.5;
    //Compute cdfvalue
    mu = size1 * size2 / (double)2.0;
    sigma = sqrt(size1 * size2 * (totalnum + 1) / (double)12.0);
    cdfvalue = normcdf(U, mu, sigma);
    delete []prank;
    delete []pflag;
    delete []ptiesize;
    return cdfvalue;
}

void getpvalue(double *pmatrix, int width, int height, std::vector<cdcps> &proregion, int size, double *pvalue)
{
    int i, j, m, n, start, end, diaelenum, trielenum, lower, upper;
    double avevalue, sumvalue, standev;
    double *ptmpmatrix, *pdiamatrix, *pupdntriangle;

    ptmpmatrix = new (std::nothrow) double[width * height];
    if(ptmpmatrix == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    else
        memcpy(ptmpmatrix, pmatrix, sizeof(double)*width*height);
    pdiamatrix = new (std::nothrow) double[size * size];
    if(pdiamatrix == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    pupdntriangle = new (std::nothrow) double[size * size];
    if(pupdntriangle == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < 2 * size; i++)
    {
        sumvalue = 0;
        for(j = i + 1; j < width; j++)
            sumvalue += *(ptmpmatrix + (j - i - 1) * width + j);
        avevalue = sumvalue / (width - i - 1);
        sumvalue = 0;
        for(j = i + 1; j < width; j++)
            sumvalue += pow(*(ptmpmatrix + (j - i - 1) * width + j) - avevalue, 2);
        standev = sqrt(sumvalue / (width - i - 2));
        for(j = i + 1; j < width; j++)
            *(ptmpmatrix + (j - i - 1) * width + j) = (*(ptmpmatrix + (j - i - 1) * width + j) - avevalue) / standev;
    }
    for(i = 0; i < (int)proregion.size(); i++)
    {
        start = proregion.at(i).start;
        end = proregion.at(i).end;
        for(j = start; j < end; j++)
        {
            diaelenum = size * size;
            for(m = 0; m < diaelenum; m++)
                *(pdiamatrix + m) = nan("");
            for(m = 0; m < size; m++)
            {
                if(j - m >= start && j < end)
                {
                    lower = (j + 1) < end? (j + 1): end;
                    upper = (j + size) < end? (j + size): end;
                    for(n = 0; n <= upper - lower; n++)
                        *(pdiamatrix + (size - m - 1) * size + n) = *(ptmpmatrix + (j - m) * width + lower + n);
                }
            }
            trielenum = 0;
            lower = (j - size) > start? (j - size): start;
            upper = j;
            for(m = lower + 1; m <= upper; m++)
            {
                for(n = lower; n < m; n++)
                {
                    *(pupdntriangle + trielenum) =  *(ptmpmatrix + n * width + m);
                    trielenum++;
                }
            }
            lower = j + 1;
            upper = (j + size) < end? (j + size): end;
            for(m = lower + 1; m <= upper; m++)
            {
                for(n = lower; n < m; n++)
                {
                    *(pupdntriangle + trielenum) =  *(ptmpmatrix + n * width + m);
                    trielenum++;
                }
            }
            *(pvalue + j) = WilcoxRanksumTest(pdiamatrix, diaelenum, pupdntriangle, trielenum);
        }
    }
    delete []ptmpmatrix;
    delete []pdiamatrix;
    delete []pupdntriangle;
}

void bintodomain(double *plocalext, int width, int height, double *pvalue, bool bstatFilter, std::vector<cdcps> &diagcd)
{
    bool isboundary;
    int i, j;
    double predata, curdata;
    cdcps tmpdiagcd;

    if(width == height)
    {
        curdata = *plocalext;
        tmpdiagcd.type = curdata == -0.5? gap: domain;
        tmpdiagcd.start = 0;
        tmpdiagcd.stpvalue = (pvalue == NULL? 0: *pvalue);
        for(i = 0; i < height; i++)
        {
            predata = curdata;
            if(i != height -1)
                curdata = *(plocalext + i);
            else
            {
                if(*(plocalext + i) != 0.5)
                    curdata = -1;
                else
                {
                    if(predata == -0.5)
                    {
                        tmpdiagcd.end = i;
                        tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                        diagcd.push_back(tmpdiagcd);
                    }
                    else if(predata == -1)
                    {
                        tmpdiagcd.type = gap;
                        tmpdiagcd.start = i;
                        tmpdiagcd.end = i;
                        tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
                        tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                        diagcd.push_back(tmpdiagcd);
                    }
                    else
                    {
                        tmpdiagcd.end = i - 1;
                        tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i - 1));
                        diagcd.push_back(tmpdiagcd);
                        tmpdiagcd.type = gap;
                        tmpdiagcd.start = i;
                        tmpdiagcd.end = i;
                        tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
                        tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                        diagcd.push_back(tmpdiagcd);
                    }
                    return;
                }
            }
            if(curdata == predata && curdata != -1)
                continue;
            else
            {
              if(curdata == -0.5 && predata == -1)
              {
                  tmpdiagcd.type = gap;
                  tmpdiagcd.start = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
              }
              else if(curdata == -0.5)
              {
                  tmpdiagcd.end = i - 1;
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i -1));
                  diagcd.push_back(tmpdiagcd);
                  tmpdiagcd.type = gap;
                  tmpdiagcd.start = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
              }
              else if((curdata == 0 || curdata == 1) && predata == -1)
              {
                  tmpdiagcd.type = domain;
                  tmpdiagcd.start = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
              }
              else if((curdata == 0 || curdata == 1) && predata == -0.5)
              {
                  tmpdiagcd.end = i - 1;
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i - 1));
                  diagcd.push_back(tmpdiagcd);
                  tmpdiagcd.type = domain;
                  tmpdiagcd.start = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
              }
              else if(curdata == 0 || curdata == 1)
                  continue;
              else if(curdata == -1 && predata == -1)
              {
                  tmpdiagcd.type = domain;
                  tmpdiagcd.start = i;
                  tmpdiagcd.end = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                  diagcd.push_back(tmpdiagcd);
              }
              else if(curdata == -1 && predata == -0.5)
              {
                  tmpdiagcd.end = i - 1;
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i - 1));
                  diagcd.push_back(tmpdiagcd);
                  tmpdiagcd.type = domain;
                  tmpdiagcd.start = i;
                  tmpdiagcd.end = i;
                  tmpdiagcd.stpvalue = (pvalue == NULL? 0: *(pvalue + i));
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                  diagcd.push_back(tmpdiagcd);
              }
              else if(curdata == -1)
              {
                  tmpdiagcd.end = i;
                  tmpdiagcd.edpvalue = (pvalue == NULL? 0: *(pvalue + i));
                  diagcd.push_back(tmpdiagcd);
              }
            }
        }
        if(bstatFilter == true)
        {
            for(i = 0; i < (int)diagcd.size(); i++)
            {
                if(diagcd.at(i).type == domain)
                {
                    isboundary = true;
                    for(j = diagcd.at(i).start; j <= diagcd.at(i).end + 1; j++)
                    {
                        if(j >= height)
                            break;
                        if(*(pvalue + j) < 0.05)
                            continue;
                        else
                        {
                            isboundary = false;
                            break;
                        }
                    }
                    if(isboundary == true)
                        diagcd.at(i).type = boundary;
                }
            }
        }
        /*//Adjuest topological domains computed in C++, so that the result can be completely the same as that of the orignal TopDom
        for(i = 0; i < (int)diagcd.size(); i++)
        {
            diagcd.at(i).start++;
            if(diagcd.at(i).end < height - 1)
            {
                diagcd.at(i).end = diagcd.at(i + 1).start + 1;
                diagcd.at(i).edpvalue = diagcd.at(i + 1).stpvalue;
            }
            else
            {
                diagcd.at(i).end = height;
                diagcd.at(i).edpvalue = diagcd.at(i).edpvalue;
            }
        }*/
    }
}

void getdiagchangepoint(double *pmatrix, int width, int height, int winsize, bool bstatFilter, std::vector<cdcps> &diagcd)
{
    int i, minsize;
    double *pdiamean, *plocalext, *pvalue;
    std::vector<int> gapidx;
    std::vector<cdcps> proregion;

    minsize = 3;
    if(width == height)
    {
        pdiamean = new (std::nothrow) double[height];
        if(pdiamean == NULL)
        {
            Errlog timerr;
            timerr.errlog("Application for memory failure;");
            exit(EXIT_FAILURE);
        }
        plocalext = new (std::nothrow) double[height];
        if(plocalext == NULL)
        {
            Errlog timerr;
            timerr.errlog("Application for memory failure;");
            exit(EXIT_FAILURE);
        }
        else
        {
            for(i = 0; i < height; i++)
                *(plocalext + i) = -0.5;
        }
        if(bstatFilter == true)
        {
            pvalue = new (std::nothrow) double[height];
            if(pvalue == NULL)
            {
                Errlog timerr;
                timerr.errlog("Application for memory failure;");
                exit(EXIT_FAILURE);
            }
            else
            {
                for(i = 0; i < height; i++)
                    *(pvalue + i) = 1;
            }
        }
        else
            pvalue = NULL;
        getdiamondmean(pmatrix, width, height, winsize, pdiamean);
        getgapindex(pmatrix, width, height, gapidx);
        getproregion(gapidx, width, height, minsize, proregion);
        detectextreme(pdiamean, width, proregion, plocalext);
        if(bstatFilter == true)
        {
            getpvalue(pmatrix, width, height, proregion, winsize, pvalue);
            for(i = 0; i < height; i++)
            {
                if(*(plocalext + i) == -1)
                {
                    if(*(pvalue + i) >= 0.05)
                        *(plocalext + i) = 0;
                }
            }
        }
        bintodomain(plocalext, width, height, pvalue, bstatFilter, diagcd);
        delete []pdiamean;
        delete []plocalext;
        delete []pvalue;
    }  
}


