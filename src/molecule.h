#pragma once

class MOLECULE
{
private:
    double **x_internal;
    int ibegin, iend;
    double innertia1,innertia2,innertia3;
    // void getRotMX();
    double dt,hdt,qdt;
    double **tmpmx1,**tmpmx2,**tmpmx3,**tmpmx4,**tmpmx5;
    void rotateInternal(); //rotaion x_internal about y for 90 deg, useful to avoid singularities.
public:
    int nat;
    double com[3],r_ms;
    double phi,theta,psi;
    double fcom[3],vcom[3];
    double torq[3],angularl[3],angularL[3]; //torque, angular momentum in space cord; angular momentum in self cord
    double **rotmx;
    void init(int na, double **xi, int i_o, double *innert, int i1, int i2);
    inline void writex();
    void computeEuler(double *xo1, double *xo2, double *xi1, double *xi2);
    void computeMolForce();
    void setDt(double ddt);
    void updateX();
    void updateV();
    void nomalizeRotmx();
    void scaleV(double sc);
    double kineticErg();
    void computeRotMK(double **tmx, double **mmx, double **smx);// tmx (dx/dq) should be 3*(3*nat), mmx (mass matrix) and smx (d2x/dq2 . F) should be 3*3
    MOLECULE();
    ~MOLECULE();
};
