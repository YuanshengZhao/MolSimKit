#include "lj.h"
#include "memory.h"
#include "global.h"
#include <cmath>

LJ::LJ(int ntp,double cf) :
fourEpsilon(nullptr), sigma(nullptr), sigmasq(nullptr)
{
    cutoff=cf;
    ntyp=ntp;
    create2DArray(fourEpsilon,ntyp,ntyp);
    create2DArray(sigma,ntyp,ntyp);
    create2DArray(sigmasq,ntyp,ntyp);
    create2DArray(shift,ntyp,ntyp);
    for(int i=0;i<ntyp;++i)
        for(int j=0;j<ntyp;++j)
            shift[i][j]=shift[j][i]=fourEpsilon[i][j]=sigmasq[i][j]=0;
}

LJ::~LJ()
{
    destroy2DArray(fourEpsilon);
    destroy2DArray(sigma);
    destroy2DArray(sigmasq);
    destroy2DArray(shift);
}

void LJ::setParam(int it,int jt,double ep,double sg)
{
    if(it<0 || it>=ntyp || jt<0 ||jt>=ntyp) throw std::runtime_error("LJ::setParam -> invalid type");
    fourEpsilon[it][jt]=fourEpsilon[jt][it]=ep*4;
    sigma[it][jt]=sigma[jt][it]=sg;
    sigmasq[it][jt]=sigmasq[jt][it]=sg*sg;
    double sr6=cub(sqr(sg/cutoff));
    shift[it][jt]=shift[jt][it]=sr6*(1-sr6);
}

void LJ::mixParam(char ty)
{
    double sr6;
    for(int i=0;i<ntyp;++i)
        for(int j=i+1;j<ntyp;++j)
        {
            fourEpsilon[i][j]=fourEpsilon[j][i]=std::sqrt(fourEpsilon[i][i]*fourEpsilon[j][j]);
            if (ty=='G')
            {
                sigma[i][j]=sigma[j][i]=std::sqrt(sigma[i][i]*sigma[j][j]);
                sigmasq[i][j]=sigmasq[j][i]=std::sqrt(sigmasq[i][i]*sigmasq[j][j]);
                sr6=cub(sqr(sigma[i][j]/cutoff));
                shift[i][j]=shift[j][i]=sr6*(1-sr6);

            }
            else if (ty=='A')
            {
                sigma[i][j]=sigma[j][i]=0.5*(sigma[i][i]+sigma[j][j]);
                sigmasq[i][j]=sigmasq[j][i]=sqr(sigma[i][j]);
                sr6=cub(sqr(sigma[i][j]/cutoff));
                shift[i][j]=shift[j][i]=sr6*(1-sr6);
            }
            else throw std::runtime_error("LJ::mixParam -> invalid mix type");
        }
}

double LJ::computeShort(int i,int j,double rsq, double *erg)
{
    double sdv6=cub(sigmasq[typ[i]][typ[j]]/rsq),spl=special[i][j];
    if(erg) *erg+=fourEpsilon[typ[i]][typ[j]]*(sdv6*(sdv6-1)+shift[typ[i]][typ[j]])*spl;
    return fourEpsilon[typ[i]][typ[j]]*(12.0/rsq)*(sdv6*(sdv6-.5))*spl;
}

double LJ::computeShortPotential(int i,int j,double rsq)
{
    double sdv6=cub(sigmasq[typ[i]][typ[j]]/rsq),spl=special[i][j];
    return fourEpsilon[typ[i]][typ[j]]*(sdv6*(sdv6-1)+shift[typ[i]][typ[j]])*spl;
}

double LJ::computeShortDuDs(int i,int j,double rsq, double *ds)
{
    double sdv6=cub(sigmasq[typ[i]][typ[j]]/rsq),spl=special[i][j];
    *ds+=fourEpsilon[typ[i]][typ[j]]*6.0/rsq*(sdv6*(.5-sdv6))*spl;
    return fourEpsilon[typ[i]][typ[j]]*42.0/sqr(rsq)*(sdv6*(sdv6-(12.0/42.0)))*spl;
}

