#include "linalg.h"
#include "../lib/OpenBLAS/include/lapacke.h"
#include "../lib/OpenBLAS/include/cblas.h"

int lin_inverse(double **mx, int n)
{
    int *ipiv=new int[n], res;
    res =LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n,n,mx[0],n,ipiv);
    res|=LAPACKE_dgetri(LAPACK_ROW_MAJOR,n,mx[0],n,ipiv);
    delete[] ipiv;
    return res;
}
int lin_inverse_sym_pos(double **mx, int n)
{
    int res=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U',n,mx[0],n);
    res|=LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U',n,mx[0],n);
    return res;
}

int lin_eigenvalue_sym(double **mx, int n, double *eigvals)
{
    double *lp_e=new double[n], *lp_tau=new double[n];
    int res=LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',n,mx[0],n,eigvals,lp_e,lp_tau);
    res|=LAPACKE_dsterf(n,eigvals,lp_e);
    delete[] lp_e;
    delete[] lp_tau;
    return res;
}

int lin_eigensystem_sym(double **mx, int n, double *eigvals)
{
    double *lp_e=new double[n], *lp_tau=new double[n];
    int res=LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',n,mx[0],n,eigvals,lp_e,lp_tau);
    res|=LAPACKE_dorgtr(LAPACK_ROW_MAJOR,'U',n,mx[0],n,lp_tau);
    res|=LAPACKE_dsteqr(LAPACK_ROW_MAJOR,'V',n,eigvals,lp_e,mx[0],n);
    delete[] lp_e;
    delete[] lp_tau;
    return res;
}

int lin_eigensystem_sym(double **mx, int n, double *eigvals, double **eigvecs)
{
    double *lp_e=new double[n], *lp_tau=new double[n];
    int res=LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',n,mx[0],n,eigvals,lp_e,lp_tau);
    res|=LAPACKE_dstevd(LAPACK_ROW_MAJOR,'V',n,eigvals,lp_e,eigvecs[0],n);
    res|=LAPACKE_dormtr(LAPACK_ROW_MAJOR,'L','U','N',n,n,mx[0],n,lp_tau,eigvecs[0],n);
    delete[] lp_e;
    delete[] lp_tau;
    return res;
}

// void lin_my_mxmx(double **mx1, double **mx2, int m,int n,int k, double **mx_o) //A_mk B_kn
// {
//     double *mmx1,mmxo;
//     for(int ii=0;ii<m;++ii)
//     {
//         mmx1=mx1[ii];
//         for(int jj=0;jj<n;++jj)
//         {
//             mmxo=0;
//             for(int kk=0;kk<k;++kk) mmxo+=mmx1[kk]*mx1[kk][jj];
//             mx_o[ii][jj]=mmxo;
//         }
//     }
// }