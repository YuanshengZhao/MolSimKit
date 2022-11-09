#pragma once

class SQDIRECT
{
private:
    int *vend,*vbegin, totalq;
    double **gvec,*qs, *gvec_first; // track the first point because the pointers in gvec does not preserve the order!!!
    double precesion;
    double dq,q_max;
    //quick sort
    void sortQ(int lo, int hi);
    //partition for quick sort
    int partitionQ(int lo, int hi);
public:
    int nbin;
    double **sq,**dsq;
    SQDIRECT(double qm, int nq, double prec);
    ~SQDIRECT();
    // compute | sum_i f_i exp(-i Q . r_i) |^2
    void computeWScalar(double **fc, int nn);
    // compute | sum_i (t_i . Q)/Q exp(-i Q . r_i) |^2, ti has shape 3*natom
    void computeWVector(double **ti, int nn);
};


