#pragma once

class PAIR_BASE
{
public:
    virtual double computeShort(int i,int j,double rsq,double *erg) = 0; // 1/r*du/dr
    virtual double computeShortPotential(int i,int j,double rsq) = 0;    // u(r)
    virtual double computeShortDuDs(int i,int j,double rsq, double *ds) = 0;         // du/ds and d2u/ds2 with s=rsq, ds2 is returned and ds is ADDED to ds
};
