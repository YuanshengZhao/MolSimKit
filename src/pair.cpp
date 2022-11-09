#include "pair.h"
#include "global.h"
#include "neighlist.h"
#include <cstdarg>
#include <cstring>
#include <iostream>

PAIR::PAIR(double cut, NEIGHLIST *ls) :
pairs(nullptr), n_pairs(0), list(nullptr)
{
    cutsq=cut*cut;
    list=ls;
}

PAIR::~PAIR()
{
    if(n_pairs) delete[] pairs;
}

void PAIR::registerPair(int nf, ...)
{
    if(n_pairs) 
    {
        std::cerr<<"PAIR::registerPair -> "<<n_pairs<<" original pairs are removed!\n";
        delete[] pairs;
    }
    n_pairs=nf;
    pairs = new PAIR_BASE*[n_pairs];
    std::va_list args;
    va_start(args,nf);
    for (int i=0; i<n_pairs; ++i)
        pairs[i]=va_arg(args,PAIR_BASE*);
    va_end(args);
}

double PAIR::compute()
{
    throw std::runtime_error("PAIR::compute -> Incomplete implementation");
    double erg=0,fpair,fp[3];
    int inum, jj, *ilist;
    double *xi,*xj,*fi,*fj,*gx,**igx,dx[3],rsq;
    for(int ii=0;ii<natom; ++ii)
    {
        inum=list->num_neigh[ii];
        ilist=list->nei_list[ii];
        igx=list->nei_dx[ii];
        xi=x[ii];
        fi=f[ii];
        for(int j=0;j<inum;++j)
        {
            jj=ilist[j];
            xj=x[jj];
            fj=f[jj];
            gx=igx[j];
            dx[0]=xi[0]-xj[0]-gx[0];
            dx[1]=xi[1]-xj[1]-gx[1];
            dx[2]=xi[2]-xj[2]-gx[2];
            if((rsq=sqr(dx[0])+sqr(dx[1])+sqr(dx[2]))<cutsq)
            {
                fpair=0;
                for(int fn=0;fn<n_pairs;++fn)
                    fpair+=pairs[fn]->computeShort(ii,jj,rsq,&erg);
                fp[0]=fpair*dx[0]; fp[1]=fpair*dx[1]; fp[2]=fpair*dx[2]; 
                fi[0]+=fp[0]; fi[1]+=fp[1]; fi[2]+=fp[2]; 
                fj[0]-=fp[0]; fj[1]-=fp[1]; fj[2]-=fp[2]; 
            }
        }
    }
    return erg;
}

double PAIR::computePotential()
{
    double erg=0;
    int inum, jj, *ilist;
    double *xi,*xj,*gx,**igx,dx[3],rsq;
    for(int ii=0;ii<natom; ++ii)
    {
        inum=list->num_neigh[ii];
        ilist=list->nei_list[ii];
        igx=list->nei_dx[ii];
        xi=x[ii];
        for(int j=0;j<inum;++j)
        {
            // std::cerr<<inum<<std::endl;
            jj=ilist[j];
            xj=x[jj];
            gx=igx[j];
            dx[0]=xi[0]-xj[0]-gx[0];
            dx[1]=xi[1]-xj[1]-gx[1];
            dx[2]=xi[2]-xj[2]-gx[2];
            if((rsq=sqr(dx[0])+sqr(dx[1])+sqr(dx[2]))<cutsq)
            {
                for(int fn=0;fn<n_pairs;++fn)
                    erg+=pairs[fn]->computeShortPotential(ii,jj,rsq);
            }
        }
    }
    return erg;
}
void PAIR::computePotential(double *erg)
{
    *erg=0;
    int inum, jj, *ilist;
    double *xi,*xj,*gx,**igx,dx[3],rsq,fdr;
    for(int ii=0;ii<natom; ++ii)
    {
        inum=list->num_neigh[ii];
        ilist=list->nei_list[ii];
        igx=list->nei_dx[ii];
        xi=x[ii];
        for(int j=0;j<inum;++j)
        {
            // std::cerr<<inum<<std::endl;
            jj=ilist[j];
            xj=x[jj];
            gx=igx[j];
            dx[0]=xi[0]-xj[0]-gx[0];
            dx[1]=xi[1]-xj[1]-gx[1];
            dx[2]=xi[2]-xj[2]-gx[2];
            fdr=0;
            if((rsq=sqr(dx[0])+sqr(dx[1])+sqr(dx[2]))<cutsq)
            {
                for(int fn=0;fn<n_pairs;++fn)
                    fdr+=pairs[fn]->computeShort(ii,jj,rsq,erg);
                dx[0]*=fdr;
                dx[1]*=fdr;
                dx[2]*=fdr;
                f[jj][0]-=dx[0];
                f[jj][1]-=dx[1];
                f[jj][2]-=dx[2];
                f[ii][0]+=dx[0];
                f[ii][1]+=dx[1];
                f[ii][2]+=dx[2];
            }
        }
    }
}

// assuming that a half list has been built
// we make use of the property d^2u/dr_i^2+\sum_j d^2u/dr_idr_j = 0
// thus we only need to calc d^2u/dr_idr_j for all i!=j pairs
void PAIR::computeHessianN(double h, double **matrix)
{
    double hh=h*.5,r_hsq=1/(h*h);
    int *ilist;
    int inum;
    double *xi,*xj,*gx,**igx,dx[3],rsq;
    int jj;
    double erg, ergdij[3], ergodij[3];
    double sd0[3],sdhhp[3],sdhhm[3];
    int ti0,ti1,ti2,tj0,tj1,tj2;
    for(int ii=0;ii<natom;++ii)
    {
        inum=list->num_neigh[ii];
        ilist=list->nei_list[ii];
        igx=list->nei_dx[ii];
        xi=x[ii];

        ti0=ii*3;ti1=ti0+1;ti2=ti1+1;


        for(int j=0;j<inum;++j)
        {
            jj=ilist[j];
            xj=x[jj];
            gx=igx[j];

            dx[0]=xi[0]-xj[0]-gx[0];
            dx[1]=xi[1]-xj[1]-gx[1];
            dx[2]=xi[2]-xj[2]-gx[2];
            sd0[0]=sqr(dx[0]); sd0[1]=sqr(dx[1]); sd0[2]=sqr(dx[2]);

            if((rsq=sd0[0]+sd0[1]+sd0[2])<cutsq)
            {
                ergdij[0]=ergdij[1]=ergdij[2]=ergodij[0]=ergodij[1]=ergodij[2]=0;

                sdhhp[0]=sqr(dx[0]+hh); sdhhp[1]=sqr(dx[1]+hh); sdhhp[2]=sqr(dx[2]+hh);
                sdhhm[0]=sqr(dx[0]-hh); sdhhm[1]=sqr(dx[1]-hh); sdhhm[2]=sqr(dx[2]-hh);
                for(int fn=0;fn<n_pairs;++fn)
                {
                    erg=2*pairs[fn]->computeShortPotential(ii,jj,rsq);
                    ergdij[0]+=erg
                               -pairs[fn]->computeShortPotential(ii,jj,sqr(dx[0]+h)+sd0[1]+sd0[2])
                               -pairs[fn]->computeShortPotential(ii,jj,sqr(dx[0]-h)+sd0[1]+sd0[2]);
                    ergdij[1]+=erg
                               -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sqr(dx[1]+h)+sd0[2])
                               -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sqr(dx[1]-h)+sd0[2]);
                    ergdij[2]+=erg
                               -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sd0[1]+sqr(dx[2]+h))
                               -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sd0[1]+sqr(dx[2]-h));

                    ergodij[0]+=pairs[fn]->computeShortPotential(ii,jj,sdhhp[0]+sdhhm[1]+sd0[2])
                                +pairs[fn]->computeShortPotential(ii,jj,sdhhm[0]+sdhhp[1]+sd0[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sdhhp[0]+sdhhp[1]+sd0[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sdhhm[0]+sdhhm[1]+sd0[2]);
                    ergodij[1]+=pairs[fn]->computeShortPotential(ii,jj,sdhhp[0]+sd0[1]+sdhhm[2])
                                +pairs[fn]->computeShortPotential(ii,jj,sdhhm[0]+sd0[1]+sdhhp[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sdhhp[0]+sd0[1]+sdhhp[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sdhhm[0]+sd0[1]+sdhhm[2]);
                    ergodij[2]+=pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sdhhm[1]+sdhhp[2])
                                +pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sdhhp[1]+sdhhm[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sdhhp[1]+sdhhp[2])
                                -pairs[fn]->computeShortPotential(ii,jj,sd0[0]+sdhhm[1]+sdhhm[2]);
                }
                ergdij[0]*=r_hsq; ergdij[1]*=r_hsq; ergdij[2]*=r_hsq; 
                ergodij[0]*=r_hsq; ergodij[1]*=r_hsq; ergodij[2]*=r_hsq; 

                tj0=jj*3; tj1=tj0+1; tj2=tj1+1;
                matrix[ti0][tj0]+=ergdij[0]; matrix[ti1][tj1]+=ergdij[1]; matrix[ti2][tj2]+=ergdij[2];
                matrix[ti0][ti0]-=ergdij[0]; matrix[ti1][ti1]-=ergdij[1]; matrix[ti2][ti2]-=ergdij[2];
                matrix[tj0][tj0]-=ergdij[0]; matrix[tj1][tj1]-=ergdij[1]; matrix[tj2][tj2]-=ergdij[2];
                
                matrix[ti0][tj1]+=ergodij[0]; matrix[ti1][tj0]+=ergodij[0];
                matrix[ti0][tj2]+=ergodij[1]; matrix[ti2][tj0]+=ergodij[1];
                matrix[ti1][tj2]+=ergodij[2]; matrix[ti2][tj1]+=ergodij[2];

                matrix[ti0][ti1]-=ergodij[0]; matrix[tj0][tj1]-=ergodij[0];
                matrix[ti0][ti2]-=ergodij[1]; matrix[tj0][tj2]-=ergodij[1];
                matrix[ti1][ti2]-=ergodij[2]; matrix[tj1][tj2]-=ergodij[2];
            }

        }
    }
}

//analytival hessian
void PAIR::computeHessian(double **matrix)
{
    int *ilist;
    int inum;
    double *xi,*xj,*gx,**igx,dx[3],rsq;
    int jj;
    double gxx,gyy,gzz,gxy,gxz,gyz;
    double cxx,cyy,czz,cxy,cxz,cyz;
    double uf1,uf2;
    int ti0,ti1,ti2,tj0,tj1,tj2;
    for(int ii=0;ii<natom;++ii)
    {
        inum=list->num_neigh[ii];
        ilist=list->nei_list[ii];
        igx=list->nei_dx[ii];
        xi=x[ii];

        ti0=ii*3;ti1=ti0+1;ti2=ti1+1;

        for(int j=0;j<inum;++j)
        {
            jj=ilist[j];
            xj=x[jj];
            gx=igx[j];

            dx[0]=xi[0]-xj[0]-gx[0];
            dx[1]=xi[1]-xj[1]-gx[1];
            dx[2]=xi[2]-xj[2]-gx[2];
            gxx=dx[0]*dx[0]; gyy=dx[1]*dx[1]; gzz=dx[2]*dx[2];
            if((rsq=gxx+gyy+gzz)<cutsq)
            {
                
                gxy=dx[0]*dx[1]; gxz=dx[0]*dx[2]; gyz=dx[1]*dx[2];
                uf1=uf2=0;
                for(int fn=0;fn<n_pairs;++fn)
                    uf2+=pairs[fn]-> computeShortDuDs(ii,jj,rsq,&uf1);
                uf2*=-4;
                uf1*=-2;

                cxx=uf1+gxx*uf2; cyy=uf1+gyy*uf2; czz=uf1+gzz*uf2;
                cxy=gxy*uf2;  cxz=gxz*uf2;  cyz=gyz*uf2;  

                tj0=jj*3; tj1=tj0+1; tj2=tj1+1;
                matrix[ti0][tj0]+=cxx; matrix[ti1][tj1]+=cyy; matrix[ti2][tj2]+=czz;
                matrix[ti0][ti0]-=cxx; matrix[ti1][ti1]-=cyy; matrix[ti2][ti2]-=czz;
                matrix[tj0][tj0]-=cxx; matrix[tj1][tj1]-=cyy; matrix[tj2][tj2]-=czz;
                
                matrix[ti0][tj1]+=cxy; matrix[ti1][tj0]+=cxy;
                matrix[ti0][tj2]+=cxz; matrix[ti2][tj0]+=cxz;
                matrix[ti1][tj2]+=cyz; matrix[ti2][tj1]+=cyz;

                matrix[ti0][ti1]-=cxy; matrix[tj0][tj1]-=cxy;
                matrix[ti0][ti2]-=cxz; matrix[tj0][tj2]-=cxz;
                matrix[ti1][ti2]-=cyz; matrix[tj1][tj2]-=cyz;
            }

        }
    }
}