#pragma once
#include "../lib/OpenBLAS/include/cblas.h"
#include "../lib/OpenBLAS/include/lapacke.h"

int lin_inverse(double **mx, int n);
int lin_inverse_sym_pos(double **mx, int n);
int lin_eigenvalue_sym(double **mx, int n, double *eigvals);

// QR eigenvalue algorithm, M=S O S^T, where O is the eigenvalues
int lin_eigensystem_sym(double **mx, int n, double *eigvals);
// Divide-and-conquer eigenvalue algorithm, M=S O S^T, where O is the eigenvalues
int lin_eigensystem_sym(double **mx, int n, double *eigvals, double **eigvecs);

// inline void lin_mxmx(double *mx1, CBLAS_TRANSPOSE t1, double *mx2, CBLAS_TRANSPOSE t2, int m,int n,int k, double beta, double *mx_o) //A_mk B_kn
// {
//     cblas_dgemm(CblasRowMajor,t1,t2,m,n,k,1,mx1,k,mx2,n,beta,mx_o,n);
// }

// inline void lin_mxvec(double *mx, CBLAS_TRANSPOSE t1, double *vec, int m,int n,double beta, double *vec_o) //A_mk B_kn
// {
//     cblas_dgemv(CblasRowMajor,t1,m,n,1,mx,n,vec,1,beta,vec_o,1);
// }

inline double lin_det3(double **mx)
{
    return mx[0][0]*(mx[1][1]*mx[2][2]-mx[1][2]*mx[2][1])
            -mx[0][1]*(mx[1][0]*mx[2][2]-mx[1][2]*mx[2][0])
            +mx[0][2]*(mx[1][0]*mx[2][1]-mx[1][1]*mx[2][0]);
}
// void lin_my_mxmx(double **mx1, double **mx2, int m,int n,int k, double **mx_o);
// int lin_eigenval(double **mx, int n);
