#ifndef DIFFTADS_H
#define DIFFTAD_H

#include "datamanager.h"

class DiffTADs
{
public:
    DiffTADs();
    ~DiffTADs();

public:
    void clear();
    void init(std::vector<contactmatrix> &cmprodata, paranorm &norm, parahierTADs &hierTADs, paraFTRLreg &FTRLreg);
    void diffanalysis(std::vector<contactdomain> &hiercd, std::vector<double> &regevar2);

public:
    std::vector<contactmatrix> m_cmprodata;
    paranorm m_paranorm;
    parahierTADs m_parahierTADs;
    paraFTRLreg m_paraFTRLreg;
    DataManager m_datamanager;
};

#endif // DIFFTAD_H
