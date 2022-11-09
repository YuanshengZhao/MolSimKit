#include "ewald.h"
#include "global.h"
#include "memory.h"
#include <cmath>
#include <iostream>

// better has erfc(r/(2^.5*sigma))<<1 -> sigma<=r_c/4
// exp(-sqr(two_pi/bl*sigma*n)/2)<<11 -> n>=4/(two_pi*sigma/bl)
EWALD::EWALD(double sg, int cf) :
g_vec(nullptr), expg(nullptr), n_gv(0)
{
    sigma=sg;
    kcutoff=cf;
    std::cerr<<"EWALD::sigma = "<<sigma<<", k_cutoff = "<<kcutoff<<std::endl;
    selfCF=-electricK/(sigma*std::sqrt(two_pi));
    r_root2sigma=1/(std::sqrt(2.0)*sigma);
    twos_rpi=2*r_root2sigma/std::sqrt(pi);
    generateG();
    computeSelfPotential();
}

EWALD::~EWALD()
{
    if(n_gv)
    {
        destroy1DArray(expg);
        destroy2DArray(g_vec);
    }
}

// void EWALD::generateG()
// {
//     if(n_gv==0)
//     {
//         create1DArray(expg ,cub(kcutoff*2+1)-1);
//         create2DArray(g_vec,cub(kcutoff*2+1)-1,3);
//     }
//     n_gv=0;
//     double longCF=1/(2*cub(bl)*epsilon_0);
//     int i2,i2j2,i2j2k2;
//     double g=two_pi/bl,gx,gy,gz;
//     double g2=g*g,gs2=-sqr(g*sigma)/2;
//     for (int i=-kcutoff; i<=kcutoff; ++i)
//     {
//         i2=i*i;
//         gx=g*i;
//         for (int j=-kcutoff; j<=kcutoff; ++j)
//         {
//             i2j2=i2+j*j;
//             gy=g*j;
//             for (int k=-kcutoff; k<=kcutoff; ++k)
//             {
//                 if (i2j2==0 && k==0) continue;
//                 i2j2k2=i2j2+k*k;
//                 gz=g*k;
//                 g_vec[n_gv][0]=gx; g_vec[n_gv][1]=gy; g_vec[n_gv][2]=gz;
//                 expg[n_gv++]=std::exp(gs2*i2j2k2)/(g2*i2j2k2)*longCF;
//             }
//         }
//     }
    // if ((cub(kcutoff*2+1)-1) != n_gv) throw std::runtime_error("EWALD::generateG -> wrong number of G vectors");
// }


//only generate half of the gs
void EWALD::generateG()
{
    if(n_gv==0)
    {
        create1DArray(expg ,(cub(kcutoff*2+1)-1)/2);
        create2DArray(g_vec,(cub(kcutoff*2+1)-1)/2,3);
    }
    n_gv=0;
    double longCF=1./(cub(bl)*epsilon_0);// as only half of the k, each weight should double
    int i2,i2j2,i2j2k2;
    double g=two_pi/bl,gx,gy,gz;
    double g2=g*g,gs2=-sqr(g*sigma)/2;
    // gx>0, no restriction on gy and gz
    for (int i=1; i<=kcutoff; ++i)
    {
        i2=i*i;
        gx=g*i;
        for (int j=-kcutoff; j<=kcutoff; ++j)
        {
            i2j2=i2+j*j;
            gy=g*j;
            for (int k=-kcutoff; k<=kcutoff; ++k)
            {
                // if (i2j2==0 && k==0) continue;
                i2j2k2=i2j2+k*k;
                gz=g*k;
                g_vec[n_gv][0]=gx; g_vec[n_gv][1]=gy; g_vec[n_gv][2]=gz;
                expg[n_gv++]=std::exp(gs2*i2j2k2)/(g2*i2j2k2)*longCF;
            }
        }
    }
    // gx=0 and gy>0, no restriction on gz
    // i2=0;
    gx=0;
    for (int j=1; j<=kcutoff; ++j)
    {
        i2j2=j*j;
        gy=g*j;
        for (int k=-kcutoff; k<=kcutoff; ++k)
        {
            // if (i2j2==0 && k==0) continue;
            i2j2k2=i2j2+k*k;
            gz=g*k;
            g_vec[n_gv][0]=gx; g_vec[n_gv][1]=gy; g_vec[n_gv][2]=gz;
            expg[n_gv++]=std::exp(gs2*i2j2k2)/(g2*i2j2k2)*longCF;
        }
    }
    // gx=0 and gy=0, gz>0
    // i2=0;
    // gx=0;
    // i2j2=0;
    gy=0;
    for (int k=1; k<=kcutoff; ++k)
    {
        // if (i2j2==0 && k==0) continue;
        i2j2k2=k*k;
        gz=g*k;
        g_vec[n_gv][0]=gx; g_vec[n_gv][1]=gy; g_vec[n_gv][2]=gz;
        expg[n_gv++]=std::exp(gs2*i2j2k2)/(g2*i2j2k2)*longCF;
    }
    if ((cub(kcutoff*2+1)-1) != n_gv*2) throw std::runtime_error("EWALD::generateG -> wrong number of G vectors");
}


double EWALD::computeLongPotential()
{
    double elong=0,ss,cc,kr,*gg,epg;
    double *si=new double[natom], *ci=new double[natom];
    for(int ig=0;ig<n_gv;++ig)
    {
        gg=g_vec[ig];
        epg=expg[ig];
        ss=cc=0;
        for (int ii=0;ii<natom;++ii)
        {
            kr=(gg[0]*x[ii][0]+gg[1]*x[ii][1]+gg[2]*x[ii][2]);
            ss+=(si[ii]=q[ii]*std::sin(kr));
            cc+=(ci[ii]=q[ii]*std::cos(kr));
        }
        elong+=(ss*ss+cc*cc)*epg;
        ss*=(epg*2);
        cc*=(epg*2);
        for (int ii=0;ii<natom;++ii)
        {
            kr=cc*si[ii]-ss*ci[ii];
            f[ii][0]+=kr*gg[0];
            f[ii][1]+=kr*gg[1];
            f[ii][2]+=kr*gg[2];
        }
    }
    delete[] si;
    delete[] ci;
    return elong;
}

void EWALD::computeHessian(double **matrix)
{
    double cc,kr,*gg,epg;
    double gxx,gyy,gzz,gxy,gxz,gyz;
    double cxx,cyy,czz,cxy,cxz,cyz;
    int ti0,ti1,ti2,tj0,tj1,tj2;
    for(int ig=0;ig<n_gv;++ig)
    {
        std::cerr<<"EWALD H "<<ig<<'/'<<n_gv<<'\r';
        gg=g_vec[ig];
        epg=expg[ig];
        gxx=gg[0]*gg[0]; gyy=gg[1]*gg[1]; gzz=gg[2]*gg[2];
        gxy=gg[0]*gg[1]; gxz=gg[0]*gg[2]; gyz=gg[1]*gg[2];
        for (int ii=0;ii<natom;++ii)
        {
            ti0=ii*3;ti1=ti0+1;ti2=ti1+1;
            for (int jj=ii+1; jj<natom; ++jj)
            {
                tj0=jj*3; tj1=tj0+1; tj2=tj1+1;

                kr=(gg[0]*(x[ii][0]-x[jj][0])+gg[1]*(x[ii][1]-x[jj][1])+gg[2]*(x[ii][2]-x[jj][2]));
                cc=2*q[ii]*q[jj]*std::cos(kr)*epg;
                cxx=gxx*cc;cyy=gyy*cc;czz=gzz*cc;cxy=gxy*cc;cxz=gxz*cc;cyz=gyz*cc;

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

double EWALD::computeSelfPotential()
{
    double q2=0;
    for (int ii=0;ii<natom;++ii) q2+=sqr(q[ii]);
    return selfErg=q2*selfCF;
}

double EWALD::computeShort(int i, int j, double drsq, double *erg)
{
    double rr=std::sqrt(drsq),rs=rr*r_root2sigma,spl=special[i][j];
    double erfcc=std::erfc(rs);
    double qij=electricK*q[i]*q[j];
    if (spl==1.0)
    {
        if(erg) *erg+=qij/rr*erfcc;
        return qij*(erfcc/(rr*drsq) + twos_rpi/drsq*std::exp(-rs*rs));
    }
    spl=1-spl;
    if(erg) *erg+=qij/rr*(erfcc-spl);
    return qij*((erfcc-spl)/(rr*drsq) + twos_rpi/drsq*std::exp(-rs*rs));
}

double EWALD::computeShortPotential(int i, int j, double drsq)
{
    drsq=std::sqrt(drsq);
    return electricK*q[i]*q[j]/drsq*(std::erfc(drsq*r_root2sigma)-(1-special[i][j]));
}

double EWALD::computeShortDuDs(int i,int j,double drsq,double *ds)
{
    double rr=std::sqrt(drsq),rs=rr*r_root2sigma;
    double erfcc=(std::erfc(rs)-(1-special[i][j]))/(rr*drsq),exps=twos_rpi*std::exp(-rs*rs)/drsq;
    double qij=electricK*q[i]*q[j];
    *ds-=0.5*qij*(erfcc + exps);
    return qij*(.75*erfcc/drsq + exps*(.75+.5*rs*rs)/drsq);
}