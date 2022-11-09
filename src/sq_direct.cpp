#include "sq_direct.h"
#include "memory.h"
#include "global.h"
#include "rnd.h"
#include <cmath>
#include "linalg.h"

#define MIN_NUM_G 10

SQDIRECT::SQDIRECT(double qm, int nq, double prec):
sq(nullptr),dsq(nullptr)
{
    q_max=qm;
    nbin=nq;
    precesion=sqr(prec*natom);
    dq=qm/nq;

    int i2,i2j2,i2j2k2;
    double g=two_pi/bl,gx,gy,gz;
    int kcutoff=q_max/g;
    totalq=cub(2*kcutoff+1)-1;
    create2DArray(gvec,totalq,3);
    gvec_first=gvec[0];
    create1DArray(qs,totalq);
    int n_gv=0, cuti2=(int)sqr(q_max/g);

    for (int i=-kcutoff; i<=kcutoff; ++i)
    {
        i2=i*i;
        gx=g*i;
        for (int j=-kcutoff; j<=kcutoff; ++j)
        {
            i2j2=i2+j*j;
            gy=g*j;
            for (int k=-kcutoff; k<=kcutoff; ++k)
            {
                if (i2j2==0 && k==0) continue;
                i2j2k2=i2j2+k*k;
                if (i2j2k2>cuti2) continue;
                gz=g*k;
                gvec[n_gv][0]=gx; gvec[n_gv][1]=gy; gvec[n_gv][2]=gz;
                qs[n_gv++]=std::sqrt((double)i2j2k2)*g;
            }
        }
    }
    // if (totalq != n_gv) throw std::runtime_error("SQDIRECT::SQDIRECT -> wrong number of G vectors");
    totalq=n_gv;
    std::cerr<<"SQDIRECT::SQDIRECT "<<totalq<<" G vecs\n";
    sortQ(0,totalq);
    // for(int i=0;i<totalq;++i) std::cerr<<i<<'\t'<<qs[i]<<'\t'<<gvec[i][0]<<'\t'<<gvec[i][1]<<'\t'<<gvec[i][2]<<std::endl;

    create1DArray(vend,nbin);
    create1DArray(vbegin,nbin);
    // create1DArray(sq,nbin);
    // create1DArray(dsq,nbin);
    vbegin[0]=i2=0;
    for(int i=1;i<nbin;++i)
    {
        gx=i*dq;
        while (i2<totalq && qs[i2]<gx) ++i2;
        vbegin[i]=vend[i-1]=i2;
    }
    vend[nbin-1]=totalq;
    for(int i=0;i<nbin;++i) std::cerr<<i<<'\t'<<i*dq<<'\t'<<(i+1)*dq<<'\t'<<vbegin[i]<<'\t'<<vend[i]<<std::endl;
    for(int i=0;i<nbin;++i) 
        if(vend[i]==vbegin[i])
            std::cerr<<"SQDIRECT::SQDIRECT Warning: Q interval #"<<i<<" has no G vecs\n"; 

    // as only q^-2 will be used afterwards, we recalculate qs
    for(int i=0;i<totalq;++i) qs[i]=1.0/(sqr(gvec[i][0])+sqr(gvec[i][1])+sqr(gvec[i][2]));
}

void SQDIRECT::sortQ(int lo, int hi)
{
    if (lo+1>=hi || lo<0)  return;
    int p=partitionQ(lo,hi); 
    sortQ(lo,p);
    sortQ(p+1,hi); 
}
int SQDIRECT::partitionQ(int lo, int hi)
{
    //use median-of-3 pivot
    int i=(lo+hi)/2;
    double temp,*tepp;
    --hi;
    if (qs[i]<qs[lo])
    {
        temp=qs[lo]; qs[lo]=qs[i]; qs[i]=temp;
        tepp=gvec[lo]; gvec[lo]=gvec[i]; gvec[i]=tepp;
        // temp=gvec[lo][0]; gvec[lo][0]=gvec[i][0]; gvec[i][0]=temp;
        // temp=gvec[lo][1]; gvec[lo][1]=gvec[i][1]; gvec[i][1]=temp;
        // temp=gvec[lo][2]; gvec[lo][2]=gvec[i][2]; gvec[i][2]=temp;
    }
    if (qs[hi]<qs[lo])
    { 
        temp=qs[lo]; qs[lo]=qs[hi]; qs[hi]=temp;
        tepp=gvec[lo]; gvec[lo]=gvec[hi]; gvec[hi]=tepp;
        // temp=gvec[lo][0]; gvec[lo][0]=gvec[hi][0]; gvec[hi][0]=temp;
        // temp=gvec[lo][1]; gvec[lo][1]=gvec[hi][1]; gvec[hi][1]=temp;
        // temp=gvec[lo][2]; gvec[lo][2]=gvec[hi][2]; gvec[hi][2]=temp;
    }
    if (qs[i]<qs[hi])
    {
        temp=qs[i]; qs[i]=qs[hi]; qs[hi]=temp;
        tepp=gvec[i]; gvec[i]=gvec[hi]; gvec[hi]=tepp;
        // temp=gvec[i][0]; gvec[i][0]=gvec[hi][0]; gvec[hi][0]=temp;
        // temp=gvec[i][1]; gvec[i][1]=gvec[hi][1]; gvec[hi][1]=temp;
        // temp=gvec[i][2]; gvec[i][2]=gvec[hi][2]; gvec[hi][2]=temp;
    }

    double pivot=qs[hi];

    i=lo-1;
    for(int j=lo;j<hi;++j)
    {
        if (qs[j]<=pivot)
        {
            ++i;
            temp=qs[i]; qs[i]=qs[j]; qs[j]=temp;
            tepp=gvec[i]; gvec[i]=gvec[j]; gvec[j]=tepp;
            // temp=gvec[i][0]; gvec[i][0]=gvec[j][0]; gvec[j][0]=temp;
            // temp=gvec[i][1]; gvec[i][1]=gvec[j][1]; gvec[j][1]=temp;
            // temp=gvec[i][2]; gvec[i][2]=gvec[j][2]; gvec[j][2]=temp;
        }
    }
    ++i;
    temp=qs[i]; qs[i]=qs[hi]; qs[hi]=temp;
    tepp=gvec[i]; gvec[i]=gvec[hi]; gvec[hi]=tepp;
    // temp=gvec[i][0]; gvec[i][0]=gvec[hi][0]; gvec[hi][0]=temp;
    // temp=gvec[i][1]; gvec[i][1]=gvec[hi][1]; gvec[hi][1]=temp;
    // temp=gvec[i][2]; gvec[i][2]=gvec[hi][2]; gvec[hi][2]=temp;
    return i;
}

SQDIRECT::~SQDIRECT()
{
    gvec[0]=gvec_first;
    destroy2DArray(gvec);
    destroy1DArray(vend);
    destroy1DArray(vbegin);
    if(sq)  destroy2DArray(sq);
    if(dsq) destroy2DArray(dsq);
    destroy1DArray(qs);
}

void SQDIRECT::computeWScalar(double **fc, int nn)
{
    double *tempp,ss,cc,kr,*si,*ci,*fci,*sqi,*dsqi;
    int vb,ve,nv,temp;
    bool flag;
    si=new double[natom];
    ci=new double[natom];
    sqi=new double[nn];
    dsqi=new double[nn];
    if(sq)  destroy2DArray(sq);
    if(dsq) destroy2DArray(dsq);
    create2DArray(sq,nn,nbin);
    create2DArray(dsq,nn,nbin);
    for(int i=0;i<nbin;++i)
    {
        nv=0;
        for(int n=0;n<nn;++n) sqi[n]=dsqi[n]=0;
        vb=vbegin[i];
        ve=vend[i];
        for(int vv=vb;vv<ve;++vv)
        {
            ++nv;
            temp=randIntStrict(ve-vv)+vv;
            tempp=gvec[vv];gvec[vv]=gvec[temp];gvec[temp]=tempp;
            cc=qs[vv];qs[vv]=qs[temp];qs[temp]=cc;
            tempp=gvec[vv];
            for(int ii=0;ii<natom;++ii)
            {
                kr=(tempp[0]*x[ii][0]+tempp[1]*x[ii][1]+tempp[2]*x[ii][2]);
                si[ii]=std::sin(kr);
                ci[ii]=std::cos(kr);
            }
            for(int n=0;n<nn;++n)
            {
                fci=fc[n];
                ss=cc=0;
                for(int ii=0;ii<natom;++ii)
                {
                    // std::cout<<fci[ii]<<'\t';
                    ss+=fci[ii]*si[ii];
                    cc+=fci[ii]*ci[ii];
                }
                // std::cout<<std::endl;
                kr=ss*ss+cc*cc;
                sqi[n]+=kr;
                dsqi[n]+=sqr(kr);
            }
            // check precision
            if(nv>=MIN_NUM_G)
            {
                flag=true;
                for(int n=0;n<nn;++n)
                    if(dsqi[n]-sqr(sqi[n])/nv>(nv-1)*nv*precesion) 
                    {
                        flag=false;
                        break;
                    }
                if(flag) break;
            }
        }
        for(int n=0;n<nn;++n)
        {
            sq[n][i]=sqi[n]/(nv*natom);
            dsq[n][i]=std::sqrt((dsqi[n]-sqr(sqi[n])/nv)/((nv-1)*nv))/natom;
        }
        std::cerr<<i<<'/'<<nbin<<','<<nv<<'\r';
    }
    std::cerr<<"SQDIRECT::computeWScalar finished\n";
    delete[] si;
    delete[] ci;
    delete[] sqi;
    delete[] dsqi;
}

void SQDIRECT::computeWVector(double **ti, int nn)
{
    double *tempp,ss,cc,kr,*si,*ci,*tii,*sqi,*dsqi,kt;
    int vb,ve,nv,temp,i_ti;
    bool flag;
    si=new double[natom];
    ci=new double[natom];
    sqi=new double[nn];
    dsqi=new double[nn];
    if(sq)  destroy2DArray(sq);
    if(dsq) destroy2DArray(dsq);
    create2DArray(sq,nn,nbin);
    create2DArray(dsq,nn,nbin);
    for(int i=0;i<nbin;++i)
    {
        nv=0;
        for(int n=0;n<nn;++n) sqi[n]=dsqi[n]=0;
        vb=vbegin[i];
        ve=vend[i];
        for(int vv=vb;vv<ve;++vv)
        {
            ++nv;
            temp=randIntStrict(ve-vv)+vv;
            tempp=gvec[vv];gvec[vv]=gvec[temp];gvec[temp]=tempp;
            cc=qs[vv];qs[vv]=qs[temp];qs[temp]=cc;
            tempp=gvec[vv];
            for(int ii=0;ii<natom;++ii)
            {
                kr=(tempp[0]*x[ii][0]+tempp[1]*x[ii][1]+tempp[2]*x[ii][2]);
                si[ii]=std::sin(kr);
                ci[ii]=std::cos(kr);
            }
            for(int n=0;n<nn;++n)
            {
                tii=ti[n];
                ss=cc=0;
                i_ti=0;
                for(int ii=0;ii<natom;++ii)
                {
                    kt=(tempp[0]*tii[i_ti]+tempp[1]*tii[i_ti+1]+tempp[2]*tii[i_ti+2]);
                    // std::cout<<kt<<'\t';
                    i_ti+=3;
                    ss+=kt*si[ii];
                    cc+=kt*ci[ii];
                }
                // std::cout<<std::endl;
                kr=(ss*ss+cc*cc)*qs[vv];
                sqi[n]+=kr;
                dsqi[n]+=sqr(kr);
            }
            // check precision
            if(nv>=MIN_NUM_G)
            {
                flag=true;
                for(int n=0;n<nn;++n)
                    if(dsqi[n]-sqr(sqi[n])/nv>(nv-1)*nv*precesion) 
                    {
                        flag=false;
                        break;
                    }
                if(flag) break;
            }
        }
        for(int n=0;n<nn;++n)
        {
            sq[n][i]=sqi[n]/(nv*natom);
            dsq[n][i]=std::sqrt((dsqi[n]-sqr(sqi[n])/nv)/((nv-1)*nv))/natom;
        }
        std::cerr<<i<<'/'<<nbin<<','<<nv<<'\r';
    }
    std::cerr<<"SQDIRECT::computeWVector finished\n";
    delete[] si;
    delete[] ci;
    delete[] sqi;
    delete[] dsqi;
}