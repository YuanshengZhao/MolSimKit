#include "neighlist.h"
#include "global.h"
#include "memory.h"
#include <cstring>

template void NEIGHLIST::build<true>(bool);
template void NEIGHLIST::build<false>(bool);

NEIGHLIST::NEIGHLIST(double cutoff, double skin, int cap) :
x_prev(nullptr), num_neigh(nullptr), nei_list(nullptr), nei_dx(nullptr), nbuild(0)
{
    create1DArray(num_neigh,natom);
    create2DArray(x_prev,natom,3);
    capacity=cap;
    create2DArray(nei_list,natom,capacity);
    create3DArray(nei_dx,natom,capacity,3);
    for(int i=0;i<natom;++i) x_prev[i][0]=x_prev[i][1]=x_prev[i][2]=1e10;
    cutsq=sqr(cutoff+skin);
    dx2=sqr(skin/2);
}

NEIGHLIST::~NEIGHLIST()
{
    destroy1DArray(num_neigh);
    destroy2DArray(x_prev);
    destroy2DArray(nei_list);
    destroy3DArray(nei_dx);
}

bool NEIGHLIST::isvalid()
{
    double *xx,*xp;
    for(int i=0;i<natom;++i)
    {
        xx=x[i];xp=x_prev[i];
        if (sqr(xx[0]-xp[0])+sqr(xx[1]-xp[1])+sqr(xx[2]-xp[2]) > dx2) return false;
    }
    return true;
}

// void NEIGHLIST::refold()
// {
//     for(int i=0;i<natom;++i)
//     {
//         x[i][0]-=std::nearbyint(x[i][0]*r_bl)*bl;
//         x[i][1]-=std::nearbyint(x[i][1]*r_bl)*bl;
//         x[i][2]-=std::nearbyint(x[i][2]*r_bl)*bl;
//     }
// }

template <bool full>
void NEIGHLIST::build(bool check)
{
    if (check && isvalid()) return;
    // refold();
    double gx[3],dx[3],mdx[3],rsq,*xx,*yy;
    int icount;
    ++nbuild;
    for(int i=0;i<natom;++i)
    {
        icount=num_neigh[i]=0;
        xx=x[i];
        for(int j = full? 0: i+1; j<natom; ++j)
        {
            if(full)
                if(i==j) continue;
            yy=x[j];
            rsq=0;
            for(int k=0;k<3;++k)
            {
                dx[k]=xx[k]-yy[k];
                gx[k]=std::nearbyint(dx[k]*r_bl)*bl;
                mdx[k]=dx[k]-gx[k];
                rsq+=sqr(mdx[k]);
            }
            if(rsq<cutsq)
            {
                nei_list[i][icount]=j;
                memcpy(nei_dx[i][icount],gx,3*sizeof(double));
                // std::cerr<<i<<"->"<<j<<' '<<icount<<std::endl;
                if((icount=++num_neigh[i])==capacity) 
                {
                    throw std::runtime_error("NEIGHLIST::build -> neighbor list overflow");
                    return;
                }
            }
        }
    }
    
    memcpy(x_prev[0],x[0],natom*(3*sizeof(double)));
}


// // cannot use parallel here!
// bool generateGhost(double cutoff)
// {
//     ntotal=natom;
//     double lim=h_bl-cutoff,nlim=-lim;
//     double *xx,dx,dy,dz;
//     bool rx[3],ry[3],rz[3];
//     rx[1]=ry[1]=rz[1]=false;
    
//     for(int i=0;i<natom;++i)
//     {
//         if (ntotal+26>=x_capacity)
//         {
//             // warning: increase typ is also needed!!!
//             // double **xn;
//             // x_capacity<<=2;
//             // create2DArray(xn,x_capacity,3);
//             // memcpy(xn[0],x[0],(x_capacity>>2)*3*sizeof(double));
//             // destroy2DArray(x);
//             // x=xn;
//             throw std::runtime_error("x_capacity too small\n");
//             return false;
//         }
//         xx=x[i];
//         rx[0]=x[i][0]<lim; rx[2]=x[i][0]>-lim;
//         ry[0]=x[i][1]<lim; ry[2]=x[i][1]>-lim;
//         rz[0]=x[i][2]<lim; rz[2]=x[i][2]>-lim;
//         for(int sx=0;sx<3;++sx)
//         {
//             if (rx[i]) continue;
//             dx=(sx-1)+bl;
//             for(int sy=0;sy<3;++sy)
//             {
//                 if(ry[i]) continue;
//                 dy=(sy-1)+bl;
//                 for(int sz=0;sz<3;++sz)
//                 {
//                     if(rz[i] || (sx==1 && sy==1 && sz==1)) continue;
//                     dz=(sz-1)+bl;
//                     x[ntotal][0]=x[i][0]+dx;
//                     x[ntotal][1]=x[i][1]+dy;
//                     x[ntotal][2]=x[i][2]+dz;
//                     g_orig[ntotal]=i;
//                     g_dx[ntotal][0]=dx;
//                     g_dx[ntotal][1]=dy;
//                     g_dx[ntotal][2]=dz;
//                     ++ntotal;
//                 }
//             }
//         }
//     }
// }

// bool updateGhost()
// {
//     double *xx,*xo,*dd;
//     for(int i=natom;i<ntotal;++i)
//     {
//         xo=x[g_orig[i]]; xx=x[i]; dd=g_dx[i];
//         xx[0]=xo[0]+dd[0];
//         xx[1]=xo[1]+dd[1];
//         xx[2]=xo[2]+dd[2];
//     }
// }