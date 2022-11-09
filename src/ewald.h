#pragma once
#include "pair_base.h"

class EWALD : public PAIR_BASE
{
private:
    double sigma;
    int kcutoff;
    double selfCF;
    double r_root2sigma,twos_rpi;
    double **g_vec;
    double *expg;
    int n_gv;
public:
    double selfErg;
    EWALD(double sg, int cf);
    ~EWALD();
    void generateG();
    double computeLongPotential();
    double computeSelfPotential();
    virtual double computeShort(int i, int j, double drsq, double *erg) override;
    virtual double computeShortPotential(int i, int j, double drsq) override;
    void computeHessian(double **matrix); // the matrix must be created with create2DArray with shape natom x natoms
    virtual double computeShortDuDs(int i,int j,double rsq,double *ds) override;
};

