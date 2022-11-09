#pragma once
#include <cmath>
#include "molecule.h"
extern double **x, **v, **f;
extern double *q, *r_m;
extern double **special;
extern int *typ;
extern int natom, ntype;//, ntotal, x_capacity;
extern int nmol;
extern MOLECULE *molecule;
extern const double epsilon_0, electricK;
extern double bl, r_bl, h_bl;
// extern double **g_dx;
// extern int *g_orig;

extern const double pi;
extern const double two_pi,four_pi;
extern const double accel;

// the second line should contain box lo and hi infomation
void initialize(int ntp, const char *xyzfile, const char *elms[], const double *r_mass, const double *charges);

// must use principle axis at COM!!!
// zero COM is not needed as it will be automatically shifted
int readMolecule(const char*xyzfile, const char *elms[], const double *r_mass, double **&xin, double *innertia);
void groupMol(int na, int nm, int ibegin, double **xi, double *innertia, int i1,int i2,double dt=0);

void inline set_bl(double bl_in){bl=bl_in; r_bl=1/bl; h_bl=bl*.5;}
template <typename numtype> inline numtype sqr(numtype x) {return x*x;}
template <typename numtype> inline numtype cub(numtype x) {return x*x*x;}
double inline distsq(const double *&x1, const double *&x2)
{
    double dx,rsq=0;
    for (int i=0;i<3;++i)
    {
        dx=std::fabs(x1[0]-x2[0]);
        dx-= static_cast<int>(dx*r_bl + 0.5)*bl;
        rsq+=dx*dx;
    }
    return rsq;
}

void printMatrix(double **mx,int m,int n,const char* info, char fmt=' ');
void writeArrayBin(double *arr,int n,const char* fname);