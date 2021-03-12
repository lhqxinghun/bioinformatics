#include <Rcpp.h>
#include <numeric>
#include <stdio.h>
#include "diffTADs.h"

// [[Rcpp::export]]
Rcpp::List rundiffTADs(Rcpp::List &matlistcon1, Rcpp::List &matlistcon2, Rcpp::String namecon1, Rcpp::String namecon2, Rcpp::List parnorm, Rcpp::List parhierTADs, Rcpp::List parFTRLreg)
{
  char nameitem[MAXCHAR];
  int i, j, m, n, k, contextnum, lencon1, lencon2, itemidx;
  double rsquared, maxrsquared, paraopt[5];
  contactmatrix tmpcmdata;
  std::vector<contactmatrix> cmprodata;
  paranorm norm;
  parahierTADs hierTADs;
  paraFTRLreg FTRLreg;
  std::vector<contactdomain> hiercd;
  std::vector<double> regevar2;
  DiffTADs DiffTADsobj;
  
  contextnum = 0;
  memset(paraopt, 0, sizeof(double) * 5);
  lencon1 = matlistcon1.length();
  lencon2 = matlistcon2.length();
  if(lencon1 > 0)
  {
    contextnum++;
    for(itemidx = 0; itemidx < lencon1; itemidx++)
    {
      Rcpp::NumericMatrix tmpmat  = matlistcon1[itemidx];
      strcpy(tmpcmdata.context, namecon1.get_cstring());
      tmpcmdata.type = rectangle;
      tmpcmdata.height = tmpmat.nrow();
      tmpcmdata.width = tmpmat.ncol();
      tmpcmdata.pvalue = (double *) new (std::nothrow) double[tmpcmdata.width * tmpcmdata.height];
      for(i = 0; i < tmpcmdata.height; i++)
      {
        for(j = 0; j < tmpcmdata.width; j++)
        {
          *(tmpcmdata.pvalue + i * tmpcmdata.width + j) = tmpmat(i,j);
        }
      }
      cmprodata.push_back(tmpcmdata);
    }
  }
  if(lencon2 > 0)
  {
    contextnum++;
    for(itemidx = 0; itemidx < lencon2; itemidx++)
    {
      Rcpp::NumericMatrix tmpmat  = matlistcon2[itemidx];
      strcpy(tmpcmdata.context, namecon2.get_cstring());
      tmpcmdata.type = rectangle;
      tmpcmdata.height = tmpmat.nrow();
      tmpcmdata.width = tmpmat.ncol();
      tmpcmdata.pvalue = (double *) new (std::nothrow) double[tmpcmdata.width * tmpcmdata.height];
      for(i = 0; i < tmpcmdata.height; i++)
      {
        for(j = 0; j < tmpcmdata.width; j++)
        {
          *(tmpcmdata.pvalue + i * tmpcmdata.width + j) = tmpmat(i,j);
        }
      }
      cmprodata.push_back(tmpcmdata);
    }
  }
  if(cmprodata.size() > 0 && contextnum > 0)
  {
    //
    norm.bmedifilter = Rcpp::as<bool> (parnorm["bmedifilter"]);
    norm.bicnorm = Rcpp::as<bool> (parnorm["bicnorm"]);
    norm.maxiter = Rcpp::as<int> (parnorm["maxiter"]);
    //
    hierTADs.winsize = Rcpp::as<int> (parhierTADs["winsize"]);
    hierTADs.bstatfilter =Rcpp::as<bool> (parhierTADs["bstatfilter"]);
    hierTADs.bTADsfilter = Rcpp::as<bool> (parhierTADs["bTADsfilter"]);
    hierTADs.minsizescale = Rcpp::as<double> (parhierTADs["minsizescale"]);
    hierTADs.maxsizescale = Rcpp::as<double> (parhierTADs["maxsizescale"]);
    //
    FTRLreg.repeatnum = Rcpp::as<int> (parFTRLreg["repeatnum"]);
    FTRLreg.iteration = Rcpp::as<int> (parFTRLreg["iteration"]);
    std::vector<double> lambdavec = Rcpp::as<std::vector<double> > (parFTRLreg["lambda"]);
    std::vector<double> alphavec = Rcpp::as<std::vector<double> > (parFTRLreg["alpha"]);
    std::vector<double> betavec = Rcpp::as<std::vector<double> > (parFTRLreg["beta"]);
    std::vector<double> l1vec = Rcpp::as<std::vector<double> > (parFTRLreg["l1"]);
    std::vector<double> l2vec = Rcpp::as<std::vector<double> > (parFTRLreg["l2"]);
    //
    if(contextnum == 1)
    {
      if(alphavec.size()==1 && betavec.size()==1 && l1vec.size()==1 && l2vec.size()==1)
      {
        FTRLreg.lambda =  0;
        FTRLreg.alpha = alphavec.at(0);
        FTRLreg.beta = betavec.at(0);
        FTRLreg.l1 = l1vec.at(0);
        FTRLreg.l2 = l2vec.at(0);
        
        for(i = 0; i < (int)hiercd.size(); i++)
        {
          if(hiercd.at(i).pvalue != NULL)
            delete []hiercd.at(i).pvalue;
        }
        hiercd.clear();
        regevar2.clear();
        DiffTADsobj.clear();
        DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
        DiffTADsobj.diffanalysis(hiercd, regevar2);
      }
      else if(alphavec.size()>=1 && betavec.size()>=1 && l1vec.size()>=1 && l2vec.size()>=1)
      {
        maxrsquared = -DBL_MAX;
        for(i = 0; i < (int)alphavec.size(); i++)
        {
          for(j = 0; j < (int)betavec.size(); j++)
          {
            for(m = 0; m < (int)l1vec.size(); m++)
            {
              for(n = 0; n < (int)l2vec.size(); n++)
              {
                FTRLreg.lambda =  0;
                FTRLreg.alpha = alphavec.at(i);
                FTRLreg.beta = betavec.at(j);
                FTRLreg.l1 = l1vec.at(m);
                FTRLreg.l2 = l2vec.at(n);
                  
                for(k = 0; k < (int)hiercd.size(); k++)
                {
                  if(hiercd.at(k).pvalue != NULL)
                    delete []hiercd.at(k).pvalue;
                }
                hiercd.clear();
                regevar2.clear();
                DiffTADsobj.clear();
                DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
                DiffTADsobj.diffanalysis(hiercd, regevar2);
                rsquared = std::accumulate(regevar2.begin(), regevar2.end(), 0.0) / regevar2.size();
                if(rsquared > maxrsquared)
                {
                  maxrsquared = rsquared;
                  paraopt[0] = FTRLreg.alpha;
                  paraopt[1] = FTRLreg.beta;
                  paraopt[2] = FTRLreg.l1;
                  paraopt[3] = FTRLreg.l2;
                }
              }
            }
          }
        }
        FTRLreg.lambda =  0;
        FTRLreg.alpha = paraopt[0];
        FTRLreg.beta = paraopt[1];
        FTRLreg.l1 = paraopt[2];
        FTRLreg.l2 = paraopt[3];
        
        for(i=0; i < (int)hiercd.size(); i++)
        {
          if(hiercd.at(i).pvalue != NULL)
            delete []hiercd.at(i).pvalue;
        }
        hiercd.clear();
        regevar2.clear();
        DiffTADsobj.clear();
        DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
        DiffTADsobj.diffanalysis(hiercd, regevar2);
      }
      else
      {
        Rcpp::stop("Null Parameter.");
      }
      //
      Rcpp::NumericMatrix reshiercd(hiercd.size(), hiercd.front().conum * hiercd.front().repnum + 4);
      Rcpp::NumericVector resregevar2 =  Rcpp::NumericVector::create();
      Rcpp::NumericVector optparameters = Rcpp::NumericVector::create();
      Rcpp::CharacterVector rescdnames = Rcpp::CharacterVector::create("left", "right", "top", "bottom");
      for(i = 0; i < (int)hiercd.size(); i++)
      {
        itemidx = 0;
        reshiercd(i,itemidx++) = hiercd.at(i).left;
        reshiercd(i,itemidx++) = hiercd.at(i).right;
        reshiercd(i,itemidx++) = hiercd.at(i).top;
        reshiercd(i,itemidx++) = hiercd.at(i).bottom;
        for(j = 0; j < hiercd.front().conum; j++)
        {
          for(m = 0; m < hiercd.front().repnum; m++)
          {
            reshiercd(i,itemidx++) = *(hiercd.at(i).pvalue + j * hiercd.front().repnum + m);
          }
        }
      }
      for(i = 0; i < hiercd.front().conum; i++)
      {
        for(j = 0; j < hiercd.front().repnum; j++)
        {
          sprintf(nameitem,"B_Con%d_R%d", i+1, j+1);
          rescdnames.push_back(nameitem);
        }
      }
      Rcpp::colnames(reshiercd) = rescdnames;
      itemidx = 0;
      if(lencon1 > 0)
      {
        for(i = 0; i < lencon1; i++)
        {
          char nameitem[MAXCHAR];
          sprintf(nameitem,"R2_Con%d_R%d", 1, i+1);
          resregevar2.push_back(regevar2.at(itemidx++), nameitem);
        }
      }
      else if(lencon2 > 0)
      {
        for(i = 0; i < lencon2; i++)
        {
          char nameitem[MAXCHAR];
          sprintf(nameitem,"R2_Con%d_R%d", 2, i+1);
          resregevar2.push_back(regevar2.at(itemidx++), nameitem);
        }
      }
      optparameters.push_back(FTRLreg.alpha, "alpha");
      optparameters.push_back(FTRLreg.beta, "beta");
      optparameters.push_back(FTRLreg.l1, "l1");
      optparameters.push_back(FTRLreg.l2, "l2");
      return Rcpp::List::create(Rcpp::Named("TADs") = reshiercd, Rcpp::Named("R2") = resregevar2, Rcpp::Named("Optpar") = optparameters);
    }
    else if(contextnum == 2)
    {
      if(alphavec.size()==1 && betavec.size()==1 && l1vec.size()==1 && l2vec.size()==1)
      {
        FTRLreg.lambda = lambdavec.at(0);
        FTRLreg.alpha = alphavec.at(0);
        FTRLreg.beta = betavec.at(0);
        FTRLreg.l1 = l1vec.at(0);
        FTRLreg.l2 = l2vec.at(0);
        for(i = 0; i < (int)lambdavec.size(); i++)
        {
          FTRLreg.lambda = lambdavec.at(i);
          
          for(j = 0; j < (int)hiercd.size(); j++)
          {
            if(hiercd.at(j).pvalue != NULL)
              delete []hiercd.at(j).pvalue;
          }          
          hiercd.clear();
          regevar2.clear();
          DiffTADsobj.clear();
          DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
          DiffTADsobj.diffanalysis(hiercd, regevar2);
        }
      }
      else if(alphavec.size()>=1 && betavec.size()>=1 && l1vec.size()>=1 && l2vec.size()>=1)
      {
        maxrsquared = -DBL_MAX;
        for(i = 0; i < (int)alphavec.size(); i++)
        {
          for(j = 0; j < (int)betavec.size(); j++)
          {
            for(m = 0; m < (int)l1vec.size(); m++)
            {
              for(n = 0; n < (int)l2vec.size(); n++)
              {
                FTRLreg.lambda =  0;
                FTRLreg.alpha = alphavec.at(i);
                FTRLreg.beta = betavec.at(j);
                FTRLreg.l1 = l1vec.at(m);
                FTRLreg.l2 = l2vec.at(n);

                for(k = 0; k < (int)hiercd.size(); k++)
                {
                  if(hiercd.at(k).pvalue != NULL)
                    delete []hiercd.at(k).pvalue;
                }                
                hiercd.clear();
                regevar2.clear();
                DiffTADsobj.clear();
                DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
                DiffTADsobj.diffanalysis(hiercd, regevar2);
                rsquared = std::accumulate(regevar2.begin(), regevar2.end(), 0.0) / regevar2.size();
                if(rsquared > maxrsquared)
                {
                  maxrsquared = rsquared;
                  paraopt[0] = FTRLreg.alpha;
                  paraopt[1] = FTRLreg.beta;
                  paraopt[2] = FTRLreg.l1;
                  paraopt[3] = FTRLreg.l2;                }
              }
            }
          }
        }
        FTRLreg.alpha = paraopt[0];
        FTRLreg.beta = paraopt[1];
        FTRLreg.l1 = paraopt[2];
        FTRLreg.l2 = paraopt[3];
        for(i = 0; i < (int)lambdavec.size(); i++)
        {
          FTRLreg.lambda = lambdavec.at(i);

          for(j = 0; j < (int)hiercd.size(); j++)
          {
            if(hiercd.at(j).pvalue != NULL)
              delete []hiercd.at(j).pvalue;
          }          
          hiercd.clear();
          regevar2.clear();
          DiffTADsobj.clear();
          DiffTADsobj.init(cmprodata, norm, hierTADs, FTRLreg);
          DiffTADsobj.diffanalysis(hiercd, regevar2);
        }
      }
      else
      {
        Rcpp::stop("Null Parameter.");
      }
      //
      Rcpp::NumericMatrix reshiercd(hiercd.size(), hiercd.front().conum * hiercd.front().repnum + 4);
      Rcpp::NumericVector resregevar2 =  Rcpp::NumericVector::create();
      Rcpp::NumericVector optparameters = Rcpp::NumericVector::create();
      Rcpp::CharacterVector rescdnames = Rcpp::CharacterVector::create("left", "right", "top", "bottom");
      for(i = 0; i < (int)hiercd.size(); i++)
      {
          itemidx = 0;
          reshiercd(i,itemidx++) = hiercd.at(i).left;
          reshiercd(i,itemidx++) = hiercd.at(i).right;
          reshiercd(i,itemidx++) = hiercd.at(i).top;
          reshiercd(i,itemidx++) = hiercd.at(i).bottom;
          for(j = 0; j < hiercd.front().conum; j++)
          {
              for(m = 0; m < hiercd.front().repnum; m++)
              {
                  reshiercd(i,itemidx++) = *(hiercd.at(i).pvalue + j * hiercd.front().repnum + m);
              }
          }
      }
      for(i = 0; i < hiercd.front().conum; i++)
      {
          for(j = 0; j < hiercd.front().repnum; j++)
          {
              sprintf(nameitem,"B_Con%d_R%d", i+1, j+1);
              rescdnames.push_back(nameitem);
          }
      }
      Rcpp::colnames(reshiercd) = rescdnames;
      itemidx = 0;
      for(i = 0; i < lencon1; i++)
      {
          char nameitem[MAXCHAR];
          sprintf(nameitem,"R2_Con%d_R%d", 1, i+1);
          resregevar2.push_back(regevar2.at(itemidx++), nameitem);
      }
      for(i = 0; i < lencon2; i++)
      {
          char nameitem[MAXCHAR];
          sprintf(nameitem,"R2_Con%d_R%d", 2, i+1);
          resregevar2.push_back(regevar2.at(itemidx++), nameitem);
      }
      optparameters.push_back(FTRLreg.alpha, "alpha");
      optparameters.push_back(FTRLreg.beta, "beta");
      optparameters.push_back(FTRLreg.l1, "l1");
      optparameters.push_back(FTRLreg.l2, "l2");
      return Rcpp::List::create(Rcpp::Named("TADs") = reshiercd, Rcpp::Named("R2") = resregevar2, Rcpp::Named("Optpar") = optparameters);
    }
    else
    {
      Rcpp::stop("More than two contexts.");
    }
  }
  else
  {
    Rcpp::stop("Empty Input data.");
  }
  return NULL;
}
