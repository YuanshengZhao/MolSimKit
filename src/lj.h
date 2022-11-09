#pragma once
#include "pair_base.h"

class LJ : public PAIR_BASE
{
private:
    double **fourEpsilon;
    double **sigma;
    double **sigmasq;
    double **shift;
    double cutoff;
    int ntyp;
public:
    LJ(int ntp,double cf);
    ~LJ();
    void setParam(int it,int jt,double ep,double sg);
    void mixParam(char ty);
    virtual double computeShort(int i,int j,double rsq,double *erg) override;
    virtual double computeShortPotential(int i,int j,double rsq) override;
    virtual double computeShortDuDs(int i,int j,double rsq,double *ds) override;
};

