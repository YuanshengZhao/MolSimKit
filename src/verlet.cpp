#include "verlet.h"
#include "global.h"
#include "memory.h"
#include <cstring>

VERLET::VERLET(double ddt) :
f_prev(nullptr)
{
    dt=ddt;
    hdt=dt*.5;
    create2DArray(f_prev,natom,3);
}

void VERLET::setDt(double ddt)
{
    dt=ddt;
    hdt=dt*.5;
}

VERLET::~VERLET()
{
    destroy2DArray(f_prev);
}

//warning: need to conver f to a before calling!
void VERLET::updateX()
{
    double *xx,*vv,*ff;
    for(int i=0;i<natom;++i)
    {
        xx=x[i];
        vv=v[i];
        ff=f[i];
        xx[0]+=(vv[0]+ff[0]*hdt)*dt;
        xx[1]+=(vv[1]+ff[1]*hdt)*dt;
        xx[2]+=(vv[2]+ff[2]*hdt)*dt;
    }
    memcpy(f_prev[0],f[0],natom*3*sizeof(double));
}

//warning: need to conver f to a before calling!
void VERLET::updateV()
{
    double *vv,*ff,*fp;
    for(int i=0;i<natom;++i)
    {
        vv=v[i];
        ff=f[i];
        fp=f_prev[i];
        vv[0]+=(fp[0]+ff[0])*hdt;
        vv[1]+=(fp[1]+ff[1])*hdt;
        vv[2]+=(fp[2]+ff[2])*hdt;
    }   
}
void VERLET::f2a()
{
    for(int i=0;i<natom;++i)
    {
        f[i][0]*=r_m[i]*accel;
        f[i][1]*=r_m[i]*accel;
        f[i][2]*=r_m[i]*accel;
    }   
}
void VERLET::scaleV(double sc)
{
    for(int i=0;i<natom;++i)
    {
        v[i][0]*=sc;
        v[i][1]*=sc;
        v[i][2]*=sc;
    }   
}
