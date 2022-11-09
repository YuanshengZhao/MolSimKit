#pragma once
#include "neighlist.h"
#include "pair_base.h"

// each pair evaluate pair force divided by dist
// arg is i_type, j_type, rsq, erg_pointer
// typedef double (*PAIRFUNCTYPE)(int,int,double,double*);
// typedef double (*PEFUNCTYPE)(int,int,double);
class PAIR
{
private:
    PAIR_BASE **pairs;
    int n_pairs;
    NEIGHLIST *list;
    double cutsq;
public:
    PAIR(double cut, NEIGHLIST *ls);
    ~PAIR();
    void registerPair(int nf, ...);
    double compute();
    double computePotential();
    void computePotential(double *erg);
    void computeHessianN(double h, double **matrix); // the matrix must be created with create2DArray with shape natom x natoms
    void computeHessian(double **matrix); // the matrix must be created with create2DArray with shape natom x natoms
};

