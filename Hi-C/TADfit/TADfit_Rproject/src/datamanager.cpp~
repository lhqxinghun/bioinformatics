#include "datamanager.h"
#include "FTRL/ftrl.h"
#include "errlog.h"
#include "ICnorm/icnorm.h"
#include "TopDom/TopDom.h"
#include <limits.h>
#include <float.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include "time.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

DataManager::DataManager()
{

}

DataManager::~DataManager()
{

}

inline void DataManager::linetovector(std::string &line, std::vector<std::string> &fields, char delimiter)
{
    std::string field;
    std::string::size_type indexb,indexe;

    indexb=0;
    indexe = line.find(delimiter, indexb);
    while(indexe != std::string::npos)
    {
        field = line.substr(indexb,indexe-indexb);
        fields.push_back(field);
        indexb=indexe+1;
        indexe = line.find(delimiter, indexb);
    }
    if(indexe == std::string::npos)
    {
        field = line.substr(indexb);
        fields.push_back(field);
    }
}

int DataManager::readtable(double *&pdestdata, int &width, int &height, std::string sourfilepath, bool bnewmem)
{
    unsigned int i, row, col, allocnum;
    std::string line;
    std::vector<std::string> fields;
    double* pbuffer = NULL;

    std::ifstream fin;
    fin.open(sourfilepath, std::ifstream::in);
    if(fin.is_open())
    {
        row = 0;
        col = UINT_MAX;
        if(bnewmem == true)
        {
            allocnum = 0;
            if((pbuffer = (double *)malloc(g_cmemory * sizeof(double))) == NULL)
            {
                Errlog timerr;
                timerr.errlog("Err: Application for memory failure;");
                exit(EXIT_FAILURE);
            }
            else
            {
                pdestdata = pbuffer;
                allocnum += g_cmemory;
            }
        }
    }
    else
    {
        Errlog timerr;
        timerr.errlog("Err: File open error;");
        return EXIT_FAILURE;
    }
    while(getline(fin, line))
    {
        linetovector(line, fields, ' ');
        if(col == fields.size()|| col == UINT_MAX)
        {
            row++;
            col = (unsigned int)fields.size();
            if(bnewmem == true)
            {
                if(row * col > allocnum)
                {
                    if((pbuffer = (double *)realloc(pbuffer, (allocnum + g_cmemstep) * sizeof(double))) ==  NULL)
                    {
                        free(pdestdata);
                        pdestdata = NULL;
                        Errlog timerr;
                        timerr.errlog("Err: Additional memory failure;");
                        exit(EXIT_FAILURE);
                    }
                    else
                    {
                        pdestdata = pbuffer;
                        allocnum += g_cmemstep;
                    }
                }
            }
            for(i = 0; i < col; i++)
            {
                pdestdata[(row - 1) * col + i] = std::stod(fields.at(i));
            }
            fields.clear();
        }
        else
        {
            free(pdestdata);
            pdestdata = NULL;
            Errlog timerr;
            timerr.errlog("Err: Invalid alignment of column;");
            exit(EXIT_FAILURE);
        }
    }
    fin.close();    
    width = col;
    height = row;
    return EXIT_SUCCESS;
}

void DataManager::writetable(std::string destfilepath, double *psourdata, int width, int height)
{
    int i, j;
    char delimiter, charbuffer[20];
    std::ofstream fout;
    std::string strvalue;

    delimiter = '\n';
    fout.open(destfilepath.c_str());
    if(fout.is_open())
    {
        for(i= 0; i < height; i++)
        {
            /*//TopDom.R
            sprintf(charbuffer, "%d ", i+1);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            sprintf(charbuffer, "%d ", i+1);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            sprintf(charbuffer, "%d ", i+1);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());

            sprintf(charbuffer, "%.4f", psourdata[i * width]);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            for(j = 1; j < width; j++)
            {
                sprintf(charbuffer, " %.4f", psourdata[i * width + j]);
                strvalue = charbuffer;
                fout.write(strvalue.c_str(), strvalue.length());
            }
            fout.write(&delimiter, sizeof(char));*/
            /*//IC-Finder.m
            sprintf(charbuffer, "%.4f", psourdata[i * width]);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            for(j = 1; j < width; j++)
            {
                sprintf(charbuffer, ",%.4f", psourdata[i * width + j]);
                strvalue = charbuffer;
                fout.write(strvalue.c_str(), strvalue.length());
            }
            fout.write(&delimiter, sizeof(char));*/
            /*//FTRL.CPP
            sprintf(charbuffer, "%.4f", psourdata[i * width]);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            for(j = 1; j < width; j++)
            {
                sprintf(charbuffer, " %.4f", psourdata[i * width + j]);
                strvalue = charbuffer;
                fout.write(strvalue.c_str(), strvalue.length());
            }
            fout.write(&delimiter, sizeof(char));*/
            //Norm.CPP
            sprintf(charbuffer, "%.4f", psourdata[i * width]);
            strvalue = charbuffer;
            fout.write(strvalue.c_str(), strvalue.length());
            for(j = 1; j < width; j++)
            {
                sprintf(charbuffer, " %.4f", psourdata[i * width + j]);
                strvalue = charbuffer;
                fout.write(strvalue.c_str(), strvalue.length());
            }
            fout.write(&delimiter, sizeof(char));
        }
    }
    fout.close();
}

void DataManager::cmlog2(double *psourdata, int width, int height, double *pdestdata)
{
    int i, len;
    double curvalue;

    len = width * height;
    for(i = 0; i < len; i++)
    {
        curvalue = *(psourdata + i);
        if(curvalue != 0)
            *(pdestdata + i) = log2(curvalue + 1);
        else
            *(pdestdata + i) = 0;
    }
}

void DataManager::cmlog10(double *psourdata, int width, int height, double *pdestdata)
{
    int i, len;
    double curvalue;

    len = width * height;
    for(i = 0; i < len; i++)
    {
        curvalue = *(psourdata + i);
        if(curvalue != 0)
            *(pdestdata + i) = log10(curvalue + 1);
        else
            *(pdestdata + i) = 0;
    }
}

void DataManager::cmedianfilter(double *psourdata, int width, int height, double *pdestdata)
{
    int i, j, k, m, n, ii, jj, min;
    double temp, window[9];

    memcpy(pdestdata, psourdata, sizeof(double)*width*height);
    for (i = 1; i < height - 1; i++)
    {
        for (j = 1; j < width - 1; j++)
        {
            k = 0;
            for (ii = i - 1; ii < i + 2; ++ii)
                for (jj = j - 1; jj < j + 2; ++jj)
                    window[k++] = *(psourdata + ii * width + jj);
            for (m = 0; m < 5; m++)
            {
                min = m;
                for (n = m + 1; n < 9; n++)
                    if (window[n] < window[min])
                        min = n;
                temp = window[m];
                window[m] = window[min];
                window[min] = temp;
            }
            *(pdestdata + i * width + j) = window[4];
        }
    }
}

void DataManager::cmznormalization(double *psourdata, int width, int height, double *pdestdata)
{
    int i, len;
    double sumvalue, avevalue, standev;

    len = width * height;
    sumvalue = 0;
    for(i = 0; i < len; i++)
        sumvalue += *(psourdata + i);
    avevalue = sumvalue / len;
    sumvalue = 0;
    for(i = 0; i < len; i++)
        sumvalue += pow(*(psourdata + i) - avevalue, 2);
    standev = sqrt(sumvalue / (len - 1));
    for(i = 0; i < len; i++)
        *(pdestdata + i) = (*(psourdata + i) - avevalue) / standev;
}

void DataManager::sdnormalization(std::vector<contactmatrix> &sourdataset)
{
    int idx, i, j, width, height;
    double itemnum, sumvalue, curvalue;
    double *pmatscale;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    pmatscale = new (std::nothrow) double[height];
    if(pmatscale == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    for(idx = 0; idx < (int)sourdataset.size(); idx++)
    {
        /*for(j = 0; j < (int)height; j++)
        {
            sumvalue = 0;
            for(m = j; m < (int)height; m++)
                sumvalue += *(psourdata +  m * width + m - j);
            *(pmatscale + j) = sumvalue == 0? 1: (sumvalue / (height - j));
        }*/
        for(i = 0; i < (int)height; i++)
        {
            itemnum = 0;
            sumvalue = 0;
            for(j = i; j < (int)height; j++)
            {
                curvalue = *(sourdataset.at(idx).pvalue +  j * width + j - i);

               //curvalue = log2(*(sourdataset.at(idx).pvalue +  j * width + j - i) + 2);

                //if(curvalue > 0)
                {
                    itemnum++;
                    sumvalue += curvalue;
                }
            }
            *(pmatscale + i) = sumvalue == 0? 1: (sumvalue / itemnum);
        }
        for(i = 0; i < height; i++)
        {
            for(j = 0; j < width; j++)
            {
                *(sourdataset.at(idx).pvalue + i * width + j) /= *(pmatscale + abs(i - j));
                //*(sourdataset.at(idx).pvalue + i * width + j) = pow(2, log2(*(sourdataset.at(idx).pvalue + i * width + j)+2) / *(pmatscale + abs(i - j)));
            }
        }
    }
    delete []pmatscale;
    //writetable("sdnorm.txt", sourdataset.at(i).pvalue, sourdataset.at(i).width, sourdataset.at(i).height);
}

void DataManager::cmpretrement(std::vector<contactmatrix> &sourdataset, pretrementmethod selmethod)
{
    int i;

    switch(selmethod)
    {
    case log10m:
        for(i = 0; i < (int)sourdataset.size(); i++)
            cmlog10(sourdataset.at(i).pvalue, sourdataset.at(i).width, sourdataset.at(i).height, sourdataset.at(i).pvalue);
        break;
    case medfm:
        for(i = 0; i < (int)sourdataset.size(); i++)
            cmedianfilter(sourdataset.at(i).pvalue, sourdataset.at(i).width, sourdataset.at(i).height, sourdataset.at(i).pvalue);
        break;
    case znorm:
        for(i = 0; i < (int)sourdataset.size(); i++)
            cmznormalization(sourdataset.at(i).pvalue, sourdataset.at(i).width, sourdataset.at(i).height, sourdataset.at(i).pvalue);
        break;
    case dnorm:
            sdnormalization(sourdataset);
        break;
    default:
        break;
    }
}

void DataManager::calpseudoreference(std::vector<contactmatrix> &sourdataset, contactmatrix &pserefdataset)
{
    int i, j, m, width, height;
    double curvalue, geomulvalue;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    strcpy(pserefdataset.context, sourdataset.front().context);
    pserefdataset.type = sourdataset.front().type;
    pserefdataset.width = width;
    pserefdataset.height = height;
    for(i = 0; i < height; i++)
    {
        for(j = 0; j <= i; j++)
        {
            geomulvalue = 0;
            for(m = 0; m < (int)sourdataset.size(); m++)
            {
                curvalue = *(sourdataset.at(m).pvalue + i * width + j);
                if(curvalue > 0)
                    geomulvalue += log(curvalue);
                else if(curvalue == 0)
                {
                    geomulvalue = 0;
                    break;
                }
                else
                {
                    Errlog timerr;
                    timerr.errlog("Appearance of illegal negative values in HiC matrix;");
                    exit(EXIT_FAILURE);
                }
            }
            curvalue = geomulvalue == 0? geomulvalue: exp(geomulvalue / sourdataset.size());
            *(pserefdataset.pvalue + i * width + j) = curvalue;
            *(pserefdataset.pvalue + j * width + i) = curvalue;
        }
    }
}

//Determine diagonal change points using PSSM. The change points can be sorted simultaneously.
void DataManager::diagcptobasecp(std::vector<std::vector<changepoint> *> &diagcp, int width, int height, std::vector<int> &basecp)
{
    //check input parameters
    /*std::cout << width << " " << height <<std::endl;
    for (int i = 0; i < diagcp.size(); i++)
    {
        std::cout << "This is the boundaries of the " << i << "th replicate: " << std::endl;
        for (int j = 0; j < diagcp.at(i)->size(); j++)
        {
            std::cout << diagcp.at(i) ->at(j).diagpos << "\t";
        }
        std::cout <<std::endl;
        for (int j = 0; j < diagcp.at(i)->size(); j++)
        {
            std::cout << diagcp.at(i) ->at(j).score << "\t";
        }
        std::cout <<std::endl;
    }*/
    
    //
    /*for (int i = 0; i < diagcp.at(0)->size(); i++)
    {
        basecp.push_back(diagcp.at(0) ->at(i).diagpos);
    }*/
    
    bool bvalid;
    int i, j, pos;
    double score, maxcpf, maxcps;
    double *psm;

    if(width == height)
    {
        psm = new (std::nothrow) double[width * 3];
        if(psm == NULL )
        {
            Errlog timerr;
            timerr.errlog("Application for memory failure;");
            exit(EXIT_FAILURE);
        }
         else
            memset(psm, 0, sizeof(double) * width * 3);
        for(i = 0; i < (int)diagcp.size(); i++)
            for(j = 0; j < (int)diagcp.at(i)->size(); j++)
            {
                pos = diagcp.at(i)->at(j).diagpos;
                score = diagcp.at(i)->at(j).score;
                psm[pos]++;
                psm[width + pos] += score;
            }
        for(i = 0; i < width; i++)
            if(psm[i] == (int)diagcp.size())
            {
                psm[i] = 0;
                psm[width + i] = 0;
                psm[width * 2 + i] = 1;
                if(i - 1 >= 0)
                {
                    psm[i - 1] = 0;
                    psm[width + i - 1] = 0;
                    psm[width * 2 + i - 1] = 0;
                }
                if(i + 1 < width)
                {
                    psm[i + 1] = 0;
                    psm[width + i + 1] = 0;
                    psm[width * 2 + i + 1] = 0;
                }
            }
        i = 0;
        while(i < width)
        {
            if(psm[width * 2 + i] == 1)
            {
                i++;
                continue;
            }
            else
            {
                j = i;
                bvalid = false;
                maxcpf = -DBL_MAX;
                while(psm[j] != 0)
                {
                    if(psm[j] > maxcpf)
                        maxcpf = psm[j];
                    j++;
                }
                j = i;
                maxcps = -DBL_MAX;
                while(psm[j] != 0)
                {
                    if(psm[j] == maxcpf && psm[width + j] > maxcps)
                    {
                        maxcps = psm[width + j];
                        pos = j;
                        bvalid = true;
                    }
                    j++;
                }
                if(bvalid == true)
                {
                    psm[pos] = 0;
                    psm[width + pos] = 0;
                    psm[width * 2 + pos] = 1;
                    i = j;
                }
                else
                    i++;
            }
        }
        
        //basecp.push_back(-1);
        for(i = 0; i < width; i++)
        {
            if (i == 0 && psm[width * 2 + i] == 1)
            {
                basecp.push_back(-1);
                continue;
            }
            if(psm[width * 2 + i] == 1)
                basecp.push_back(i);
        }
        delete []psm;
    }
}


void DataManager::icnormalization(std::vector<contactmatrix> &sourdataset, int maxiter)
{
    int i, j, width, height;
    double sumvalue, refvalue;
    contactmatrix pserefdataset;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    pserefdataset.pvalue = new (std::nothrow) double[width * height];
    if(pserefdataset.pvalue == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    calpseudoreference(sourdataset, pserefdataset);
    sumvalue = 0;
    for(i = 0; i < height; i++)
    {
        for(j = 0; j < width; j++)
            sumvalue += *(pserefdataset.pvalue + i * width + j);
    }
    refvalue = sumvalue / height;
    for(i = 0; i < (int)sourdataset.size(); i++)
        normalizationbyic(sourdataset.at(i).pvalue, width, height, maxiter, refvalue, sourdataset.at(i).pvalue);
    delete []pserefdataset.pvalue;
}


//winsize >= 2
void DataManager::findiagcpbyTopDom_gmean(std::vector<contactmatrix> &sourdataset, int winsize, bool bstatFilter, std::vector<changepoint> &diagcp)
{
    int i, width, height;
    changepoint tcp;
    std::vector<cdcps> diagcd;
    contactmatrix pserefdataset;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    pserefdataset.pvalue = new (std::nothrow) double[width * height];
    if(pserefdataset.pvalue == NULL)
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    calpseudoreference(sourdataset, pserefdataset);
    getdiagchangepoint(pserefdataset.pvalue, width, height, winsize, bstatFilter, diagcd);
    //Change points for the end of each topological domain by TopDom
    diagcp.push_back(changepoint{-1, 1});
    for(i = 0; i < (int)diagcd.size(); i++)
    {
        tcp.diagpos = diagcd.at(i).end;
        if(bstatFilter == true)
        {
            if(diagcd.at(i).type == gap || diagcd.at(i).end == width - 1)
                tcp.score = 1;
            else
                tcp.score = 1- diagcd.at(i).edpvalue;
        }
        else
            tcp.score = 1;
        diagcp.push_back(tcp);
    }
    //Examing change point for the end of last topological domain
    if(diagcp.at(diagcp.size()-1).diagpos != width - 1)
    {
        tcp.diagpos = width - 1;
        tcp.score = 1;
        diagcp.push_back(tcp);
    }
    delete []pserefdataset.pvalue;
}

void DataManager::findiagcpbyTopDom_pssm(std::vector<contactmatrix> &sourdataset, int winsize, bool bstatFilter, std::vector<changepoint> &diagcp)
{
	int m, i, width, height;
    changepoint tcp;
    std::vector<changepoint> *pdiagcp;
    std::vector<std::vector<changepoint>*> pdiagcp_array;
    std::vector<int> basecp;
	//find boundaries
	for(m = 0; m < sourdataset.size(); m++)
    {
        std::vector<cdcps> diagcd;
        pdiagcp = new std::vector<changepoint>;
        //std::cout <<"processing replicate "<< m <<";" << std::endl;
        width = sourdataset.at(m).width;
        height = sourdataset.at(m).height;
        getdiagchangepoint(sourdataset.at(m).pvalue, width, height, winsize, bstatFilter, diagcd);
        pdiagcp->push_back(changepoint{0, 1});
        //std::cout <<"The number of boundaris is:"<< diagcd.size() << std::endl;
        for ( i = 0; i < diagcd.size(); i++)
        {
            tcp.diagpos = diagcd.at(i).end;
            if(bstatFilter == true)
            {
                if(diagcd.at(i).type == gap || diagcd.at(i).end == width - 1)
                    tcp.score = 1;
                else
                    tcp.score = 1- diagcd.at(i).edpvalue;
            }
        else
            tcp.score = 1;
        pdiagcp->push_back(tcp);
        }
        diagcd.clear();
        //Examing change point for the end of last topological domain
        if(pdiagcp->at(pdiagcp->size()-1).diagpos != width - 1)
        {
            tcp.diagpos = width - 1;
            tcp.score = 1;
            pdiagcp->push_back(tcp);
        }
        pdiagcp_array.push_back(pdiagcp);
    }

    /*std::cout << pdiagcp_array.size() << std::endl;
    std::cout <<"TAD boundaries for each replicate are called" << std::endl;
    for(int j = 0; j < pdiagcp_array.size(); j++)
    {
        std::cout << pdiagcp_array.at(j)->size() << std::endl;   
    }*/


    diagcptobasecp(pdiagcp_array, width, height, basecp);
    for (i = 0; i < basecp.size(); i++)
    {
        diagcp.push_back(changepoint{basecp.at(i), 1});
    }
    //std::cout <<"TAD boundaries are processed with pssm, " << "total number of TAD boundaries are" << basecp.size() << std::endl;

}


void DataManager::mapcontextindex(std::vector<contactmatrix> &sourdataset, std::vector<int> &contextoindex, std::vector<std::string> &indextocontext)
{
    int i, index;
    std::map<std::string, int> contexindex;

    for(i = 0; i < (int)sourdataset.size(); i++)
        contexindex.insert(std::pair<std::string,int>(std::string(sourdataset.at(i).context), contexindex.size()));
    for(i = 0; i < (int)sourdataset.size(); i++)
    {
        index = contexindex[std::string(sourdataset.at(i).context)];
        contextoindex.push_back(index);
        while(index >= (int)indextocontext.size())
            indextocontext.push_back(std::string());
        indextocontext[index] = std::string(sourdataset.at(i).context);
    }
}

/*
void DataManager::diagcptohiercd(std::vector<contactmatrix> &sourdataset, std::vector<changepoint> &diagcp, int repeatnum, bool bTADsfilter, double minsizescale, double maxsizescale, std::vector<contactdomain> &hiercd)
{
    bool bvalid;
    int i, j, m, n, width, conum;
    double minthre, maxthre, *tempvalue;
    contactdomain hiercditem;
    std::vector<int> contoidx;
    std::vector<std::string> idxtocon;

    width = sourdataset.front().width;
    minthre = minsizescale * (width-1);
    maxthre = maxsizescale * (width-1);
    mapcontextindex(sourdataset, contoidx, idxtocon);
    conum = (int)idxtocon.size();
    for(i = 0; i < (int)diagcp.size() - 1; i++)
    {
       for(j = 0; j <= i; j++)
       {
           for(m = i; m >= 0; m--)
           {
               for(n = j; n <= i; n++)
               {
                   if(m >= j)
                       bvalid = true;
                   else if(i > n)
                       bvalid = true;
                   else
                       bvalid = false;
                   if(bvalid == true)
                   {
                       hiercditem.left = diagcp.at(j).diagpos + 1;
                       hiercditem.right = diagcp.at(n + 1).diagpos;
                       hiercditem.top = diagcp.at(m).diagpos + 1;
                       hiercditem.bottom = diagcp.at(i + 1).diagpos;
                       hiercditem.conum = conum;
                       hiercditem.repnum = repeatnum;
                       if(bTADsfilter == true)
                       {
             if((hiercditem.bottom - hiercditem.top) > maxthre ||(hiercditem.right - hiercditem.left) > maxthre || (hiercditem.bottom - hiercditem.top) < minthre || (hiercditem.right - hiercditem.left) < minthre)
             continue;
                       }
                       tempvalue = new (std::nothrow) double[conum * repeatnum];
                       if(tempvalue == NULL )
                       {
                         Errlog timerr;
                         timerr.errlog("Application for memory failure;");
                         exit(EXIT_FAILURE);
                       }
                       else
                         memset(tempvalue, 0, sizeof(double) * conum * repeatnum);
                       hiercditem.pvalue = tempvalue;
                       hiercd.push_back(hiercditem);
                   }
               }
           }
       }
    }

   std::vector<contactdomain>::iterator it= hiercd.begin();
   char filename[30]="./hiercd.txt"; 
   std::ofstream ofs(filename);
   for(; it!=hiercd.end();it++)

   {
      ofs<<(*it).left<<"\t"<<(*it).right<<"\t"<<(*it).top<<"\t"<<(*it).bottom<<"\n";
   }
      ofs.close();
}*/


// my diagcptohiercd, only consider the traditional hierarchical TADs


void DataManager::diagcptohiercd(std::vector<contactmatrix> &sourdataset, std::vector<changepoint> &diagcp, int repeatnum, bool bTADsfilter, double minsizescale, double maxsizescale, std::vector<contactdomain> &hiercd)
{
    int i, j, width, conum;
    double minthre, maxthre, *tempvalue;
    contactdomain hiercditem;
    std::vector<int> contoidx;
    std::vector<std::string> idxtocon;

    width = sourdataset.front().width;
    minthre = minsizescale * (width-1);
    maxthre = maxsizescale * (width-1);
    mapcontextindex(sourdataset, contoidx, idxtocon);
    conum = (int)idxtocon.size();
    for(i = 0; i < (int)diagcp.size() - 1; i++)
    {
       for(j = 0; j <= i; j++)
       {

          hiercditem.left = diagcp.at(j).diagpos + 1;
          hiercditem.right = diagcp.at(i + 1).diagpos;
          hiercditem.top = diagcp.at(j).diagpos + 1;
          hiercditem.bottom = diagcp.at(i + 1).diagpos;
          hiercditem.conum = conum;
          hiercditem.repnum = repeatnum;
          if(bTADsfilter == true)
          {
             if((hiercditem.bottom - hiercditem.top ) > maxthre ||(hiercditem.right - hiercditem.left ) > maxthre || (hiercditem.bottom - hiercditem.top ) < minthre || (hiercditem.right - hiercditem.left ) < minthre)
             continue;
           }
           tempvalue = new (std::nothrow) double[conum * repeatnum];
           if(tempvalue == NULL )
           {
              Errlog timerr;
              timerr.errlog("Application for memory failure;");
              exit(EXIT_FAILURE);
           }
           else
              memset(tempvalue, 0, sizeof(double) * conum * repeatnum);
              hiercditem.pvalue = tempvalue;
              hiercd.push_back(hiercditem);             
       }
    }

   /*std::vector<contactdomain>::iterator it= hiercd.begin();
   char filename[30]="./hiercd.txt"; 
   std::ofstream ofs(filename);
   for(; it!=hiercd.end();it++)
   {
      ofs<<(*it).left<<"\t"<<(*it).right<<"\t"<<(*it).top<<"\t"<<(*it).bottom<<"\n";
   }
      ofs.close();*/
}


void DataManager::hiercdfilter(std::vector<contactmatrix> &sourdataset, double minsizescale, double maxsizescale, std::vector<contactdomain> &hiercd)
{
    bool *bfiltered;
    int i, width;
    std::vector<contactdomain> temphiercd;

    width = sourdataset.front().width;
    bfiltered = new (std::nothrow) bool[hiercd.size()];
    if(bfiltered == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    else
        memset(bfiltered, 0, sizeof(bool) * hiercd.size());
    //Size filter
    if(minsizescale > 0 || maxsizescale > 0)
    {
        for(i = 0; i < (int)hiercd.size(); i++)
        {
            if(bfiltered[i] == false)
            {
                if((hiercd.at(i).bottom - hiercd.at(i).top + 1) > width/2 || (hiercd.at(i).right - hiercd.at(i).left + 1) > width/2
                        || (hiercd.at(i).bottom - hiercd.at(i).top + 1) < 2 || (hiercd.at(i).right - hiercd.at(i).left + 1) < 2)//maxsizescale * width//minsizescale * width
                    bfiltered[i] = true;
            }
        }
    }
    //Remove filtered contact domain
    for(i = 0; i < (int)hiercd.size(); i++)
    {
        if(bfiltered[i] != true)
            temphiercd.push_back(hiercd.at(i));
    }
    hiercd.clear();
    for(i = 0; i < (int)temphiercd.size(); i++)
        hiercd.push_back(temphiercd.at(i));
    delete []bfiltered;
}

void ftrlwrite(std::vector<contactmatrix> &sourdataset, std::vector<contactdomain> &hiercd, int selcdnum, gsl_matrix *&pftrlD, gsl_matrix *pftrlY, gsl_matrix *pftrlX, gsl_matrix *pftrlB)
{
    int i, j, m, n, width, height, Ydim, Xdim, Bdim, sampidx;
    double curbvalue, curxvalue, curpvalue;
    gsl_matrix *pcalftrlB;
    gsl_matrix *pcalftrlP;
    std::ofstream fout;
    std::ostringstream strvalue;

    width = sourdataset.front().width;
    height = sourdataset.front().height;
    Ydim = (int)sourdataset.size();
    Xdim = (int)hiercd.size();
    Bdim = pftrlB->size2;
    pcalftrlB = gsl_matrix_alloc(Xdim, Bdim);
    pcalftrlP = gsl_matrix_alloc((width + 1) * height / 2, 1);
    for(i = 0; i < Xdim; i++)
    {
        for(j = 0; j < Bdim; j++)
        {
            curbvalue = fabs(gsl_matrix_get(pftrlB, i, j));
            gsl_matrix_set(pcalftrlB, i, j, curbvalue);
        }
    }
    //Write Y
    strvalue << "./out/y.out";
    fout.open(strvalue.str(), std::ofstream::out | std::ofstream::trunc); //if file exists, reset it.
    if(fout.is_open())
    {
        strvalue << std::setiosflags(std::ios::fixed) << std::setprecision(4);
        strvalue.str("");
        for(i = 0; i < Ydim; i++)
        {
            for(j = 0; j < height; j++)
            {
                for(m = 0; m <= j; m++)
                {
                    if(j == 0 && m == 0)
                        strvalue << *(sourdataset.at(i).pvalue + j * width + m);
                    else
                        strvalue << " " << *(sourdataset.at(i).pvalue + j * width + m);
                }
            }
            strvalue << std::endl;
            fout.write(strvalue.str().c_str(), strvalue.str().length());
            strvalue.str("");
        }
    }
    fout.close();
    //Write P
    strvalue << "./out/p.out";
    fout.open(strvalue.str(), std::ofstream::out | std::ofstream::app);     //if file exists, add content at the end
    if(fout.is_open())
    {
        for(i = 0; i < Bdim; i++)
        {
            sampidx = 0;
            for(j = 0; j < height; j++)
            {
                for(m = 0; m <= j; m++)
                {
                    curpvalue = 0;
                    ftrlgetsample(sourdataset, j, m, hiercd, pftrlD, pftrlY, pftrlX, selcdnum);
                    for(n = 0; n < selcdnum; n++)
                    {
                        curxvalue = gsl_matrix_get(pftrlX, 0, n);
                        if(curxvalue != 0)
                            curpvalue += gsl_matrix_get(pcalftrlB, n, i) * curxvalue;
                    }
                    gsl_matrix_set(pcalftrlP, sampidx, 0, curpvalue);
                    sampidx++;
                }
            }
            strvalue << std::setiosflags(std::ios::fixed) << std::setprecision(4);
            strvalue.str("");
            for(j = 0; j < (int)pcalftrlP->size1; j++)
            {
                if(j == 0)
                    strvalue << gsl_matrix_get(pcalftrlP, j, 0);
                else
                    strvalue << " " << gsl_matrix_get(pcalftrlP, j, 0);
            }
            strvalue << std::endl;
            fout.write(strvalue.str().c_str(), strvalue.str().length());
            strvalue.str("");
        }
    }
    fout.close();
    gsl_matrix_free(pcalftrlB);
    gsl_matrix_free(pcalftrlP);
}

//Save result of each iteration
void DataManager::ftrlsolver(std::vector<contactmatrix> &sourdataset, std::vector<contactdomain> &hiercd, std::vector<double> &regevar2,
                             int repeatnum, int iteration, double lambda, double alpha, double beta, double l1, double l2)
{
    int r, i, j, m, height, Ddim, Ydim, Xdim, Bdim;
    double curbvalue;
    std::vector<int> contoidx; //like [0,0,0,0,1,1,1,1]
    std::vector<std::string> idxtocon; //like[“GM12878”, "IMR90"]
    std::vector<std::pair<int, int>> randcoord;//random coordinate pair
    gsl_matrix *pftrlD;//
    gsl_matrix *pftrlY;//Y
    gsl_matrix *pftrlX;//X
    gsl_matrix *pftrlB;//B
    gsl_matrix *pftrlZ;//
    gsl_matrix *pftrlN;//
    gsl_matrix *pftrlP;//

    height = sourdataset.front().height;
    mapcontextindex(sourdataset, contoidx, idxtocon);
    Ddim = (int)height;
    Ydim = (int)sourdataset.size();
    Xdim = (int)hiercd.size();
    Bdim = (int)idxtocon.size();
    ftrlalloc(pftrlD, pftrlY, pftrlX, pftrlB, pftrlZ, pftrlN, pftrlP, Ddim, Ydim, Xdim, Bdim);
    for(r = 0; r < repeatnum; r++)
    {
	std::cout << "repeat num:" << r << std::endl;
        gsl_matrix_set_zero(pftrlB);//
        gsl_matrix_set_zero(pftrlZ);
        gsl_matrix_set_zero(pftrlN);
	
        for(j = 0; j < height; j++)
        {
            for(m = 0; m <= j; m++)
                randcoord.push_back(std::pair<int, int>(j, m));
        }
        random_shuffle(randcoord.begin(), randcoord.end());

	for(i = 0; i < iteration; i++)
        {
            std::cout << "iteration num:" << i << std::endl;
            for(j = 0; j < (int)randcoord.size(); j++)
            {
		
                ftrlgetsample(sourdataset, randcoord.at(j).first, randcoord.at(j).second, hiercd, pftrlD, pftrlY, pftrlX, Xdim);
                ftrlpredict(pftrlX, pftrlB, pftrlZ, pftrlN, pftrlP, alpha, beta, l1, l2);
                ftrlupdate(pftrlY, pftrlX, pftrlB, pftrlZ, pftrlN, pftrlP, contoidx, lambda, alpha);
            }
            
            ftrlwrite(sourdataset, hiercd, Xdim, pftrlD, pftrlY, pftrlX, pftrlB);

        }
	randcoord.clear();  
	
        for(i = 0; i < Xdim; i++)
        {
            for(j = 0; j < Bdim; j++)
            {
                curbvalue = fabs(gsl_matrix_get(pftrlB, i, j)); //write beta to hiercd
                *(hiercd.at(i).pvalue + j * repeatnum + r) = curbvalue;
            }
        }
    }
    ftrloutput(sourdataset, hiercd, regevar2, Xdim, pftrlD, pftrlY, pftrlX, contoidx);
    ftrlfree(pftrlD, pftrlY, pftrlX, pftrlB, pftrlZ, pftrlN, pftrlP);
}
