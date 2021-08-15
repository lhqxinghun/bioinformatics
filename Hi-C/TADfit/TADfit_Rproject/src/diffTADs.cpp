/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "diffTADs.h"
#include "errlog.h"
#include <string.h>
#include<iostream>
#include <sstream>
#include <fstream>
#include <string>
#ifdef _OPENMP
#include "omp.h"
#endif
using namespace std;

DiffTADs::DiffTADs()
{

}

DiffTADs::~DiffTADs()
{
    clear();
}

void DiffTADs::clear()
{
    int i;
    
    for(i = 0; i < (int)m_cmprodata.size(); i++)
    {
      if(m_cmprodata.at(i).pvalue != NULL)
      {
        delete []m_cmprodata.at(i).pvalue;
        m_cmprodata[i].pvalue = NULL;
      }
    }
    m_cmprodata.clear();
}

void DiffTADs::init(std::vector<contactmatrix> &cmprodata, paranorm &norm, parahierTADs &hierTADs, paraFTRLreg &FTRLreg)
{
    int i, itemnum;
    contactmatrix tmpcmdata;

    itemnum = (int)cmprodata.size();
    for(i = 0; i < itemnum; i++)
    {
        strcpy(tmpcmdata.context, cmprodata.at(i).context);
        tmpcmdata.type = (geotype)((int)cmprodata.at(i).type);
        tmpcmdata.width = cmprodata.at(i).width;
        tmpcmdata.height = cmprodata.at(i).height;
        tmpcmdata.pvalue = new (std::nothrow) double[tmpcmdata.width * tmpcmdata.height];
        if(tmpcmdata.pvalue == NULL)
        {
          Errlog timerr;
          timerr.errlog("Application for memory failure;");
          exit(EXIT_FAILURE);
        }
        memcpy(tmpcmdata.pvalue, cmprodata.at(i).pvalue, sizeof(double)*tmpcmdata.width*tmpcmdata.height);
        m_cmprodata.push_back(tmpcmdata);
    }
    memcpy(&m_paranorm, &norm, sizeof(paranorm));
    memcpy(&m_parahierTADs, &hierTADs, sizeof(parahierTADs));
    memcpy(&m_paraFTRLreg, &FTRLreg, sizeof(paraFTRLreg));
}

void DiffTADs::diffanalysis(std::vector<contactdomain> &hiercd, std::vector<double> &regevar2)
{
    bool btaskvalid;
    int i, itemnum;
    //Parameters for median filter and distance normalization
    bool bmedifilter;
    bool bicnorm;
    int maxiter;
    //Parameters for hierarchical TADs
    int winsize;
    bool bstatfilter;
    bool bTADsfilter;
    double minsizescale;
    double maxsizescale;
    //Parameters for FTRL reggression
    int repeatnum;
    int iteration;
    double lambda;
    double alpha;
    double beta;
    double l1;
    double l2;
    //
    bmedifilter = m_paranorm.bmedifilter;
    bicnorm = m_paranorm.bicnorm;
    maxiter = m_paranorm.maxiter;
    winsize = m_parahierTADs.winsize;
    bstatfilter = m_parahierTADs.bstatfilter;
    bTADsfilter = m_parahierTADs.bTADsfilter;
    minsizescale = m_parahierTADs.minsizescale;
    maxsizescale = m_parahierTADs.maxsizescale;
    repeatnum = m_paraFTRLreg.repeatnum;
    iteration = m_paraFTRLreg.iteration;
    lambda = m_paraFTRLreg.lambda;
    alpha = m_paraFTRLreg.alpha;
    beta = m_paraFTRLreg.beta;
    l1 = m_paraFTRLreg.l1;
    l2 = m_paraFTRLreg.l2;
    btaskvalid = true;
    //
    itemnum = m_cmprodata.size();
    for(i = 0; i < itemnum - 1; i++)
    {
        if((m_cmprodata.at(i).width != m_cmprodata.at(i).height) ||
                (m_cmprodata.at(i + 1).width != m_cmprodata.at(i + 1).height) || (m_cmprodata.at(i).width != m_cmprodata.at(i + 1).width))
        {
            Errlog timerr;
            timerr.errlog("Application for memory failure;");
            btaskvalid = false;
        }
    }
    //
    if(btaskvalid == true)
    {
        std::vector<changepoint> diagcp;
        if(bmedifilter == true)
            m_datamanager.cmpretrement(m_cmprodata, medfm);
        if(bicnorm == true)
            m_datamanager.icnormalization(m_cmprodata, maxiter);
        m_datamanager.findiagcpbyTopDom(m_cmprodata, winsize, bstatfilter, diagcp);

        ofstream f;
        f.open("/home/biology/liucongyi/Hic/hicda/DiffTADs/out/basetads.txt", ios::app);
	for(int i=0; i< (int)diagcp.size()-1;i++){
		for(int j=i;j>=0;j--){
		   f<<diagcp[i].diagpos+1<<"\t"<<diagcp[i+1].diagpos<<"\t"<<diagcp[j].diagpos+1<<"\t"<<diagcp[j+1].diagpos<<endl;
                    //f<<diagcp[i]+1<<endl;
		}
	}
	f.close(); 
        for(int i=0; i< (int)diagcp.size();i++){
           std::cout <<  diagcp[i].diagpos<<"  "<<diagcp[i].score<< std::endl;
        }

        m_datamanager.diagcptohiercd(m_cmprodata, diagcp, repeatnum, bTADsfilter, minsizescale, maxsizescale, hiercd);
        //if(bTADsfilter == true)
        //m_datamanager.hiercdfilter(m_cmprodata, minsizescale, maxsizescale, hiercd);
        m_datamanager.cmpretrement(m_cmprodata, log10m);
        m_datamanager.ftrlsolver(m_cmprodata, hiercd, regevar2, repeatnum, iteration, lambda, alpha, beta, l1, l2);
    }
}
