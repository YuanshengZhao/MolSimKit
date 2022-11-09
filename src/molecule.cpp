#include "molecule.h"
#include <cmath>
#include <cstring>
#include "global.h"
#include "linalg.h"
#include "memory.h"

#define ITER_ACCU 1e-10
#define ROT_THR .707107
#define MOL_AUTOROT

MOLECULE::MOLECULE():
x_internal(nullptr),dt(0),hdt(0),qdt(0),tmpmx1(nullptr),tmpmx2(nullptr),tmpmx3(nullptr),tmpmx4(nullptr),tmpmx5(nullptr)
{
	create2DArray(rotmx,3,3);
	create2DArray(tmpmx1,3,3);
	create2DArray(tmpmx2,3,3);
	create2DArray(tmpmx3,3,3);
	create2DArray(tmpmx4,3,3);
	create2DArray(tmpmx5,3,3);
}
MOLECULE::~MOLECULE()
{
	destroy2DArray(x_internal);
	destroy2DArray(rotmx);
	destroy2DArray(tmpmx1);
	destroy2DArray(tmpmx2);
	destroy2DArray(tmpmx3);
	destroy2DArray(tmpmx4);
	destroy2DArray(tmpmx5);
}

void MOLECULE::setDt(double ddt)
{
    dt=ddt;
    hdt=dt*.5;
	qdt=dt*.25;
}
void MOLECULE::scaleV(double sc)
{
	vcom[0]    *=sc;
	vcom[1]    *=sc;
	vcom[2]    *=sc;
	angularl[0]*=sc;
	angularl[1]*=sc;
	angularl[2]*=sc;	
}

void MOLECULE::nomalizeRotmx()
{
	double q1=rotmx[0][0]+rotmx[1][1]+rotmx[2][2],q2,q3,q4,temp;
	int sel=4;
	for(int i=0;i<3;++i)
		if(rotmx[i][i]>q1)
		{
			sel=i;
			q1=rotmx[i][i];
		}
	switch (sel)
	{
	case 0:
		q2=std::sqrt(1+q1-rotmx[1][1]-rotmx[2][2])*.5;
		temp=.25/q2;
		q1=(rotmx[2][1]-rotmx[1][2])*temp;
		q3=(rotmx[0][1]+rotmx[1][0])*temp;
		q4=(rotmx[2][0]+rotmx[0][2])*temp;
		break;
	case 1:
		q3=std::sqrt(1+q1-rotmx[0][0]-rotmx[2][2])*.5;
		temp=.25/q3;
		q1=(rotmx[0][2]-rotmx[2][0])*temp;
		q2=(rotmx[0][1]+rotmx[1][0])*temp;
		q4=(rotmx[1][2]+rotmx[2][1])*temp;
		break;
	case 2:
		q4=std::sqrt(1+q1-rotmx[0][0]-rotmx[1][1])*.5;
		temp=.25/q4;
		q1=(rotmx[1][0]-rotmx[0][1])*temp;
		q2=(rotmx[2][0]+rotmx[0][2])*temp;
		q3=(rotmx[2][1]+rotmx[1][2])*temp;
		break;
	default:
		q1=std::sqrt(1+q1)*.5;
		temp=.25/q1;
		q2=(rotmx[2][1]-rotmx[1][2])*temp;
		q3=(rotmx[0][2]-rotmx[2][0])*temp;
		q4=(rotmx[1][0]-rotmx[0][1])*temp;
		break;
	}
	temp=1.0/std::sqrt(sqr(q1)+sqr(q2)+sqr(q3)+sqr(q4));
	q1*=temp;
	q2*=temp;
	q3*=temp;
	q4*=temp;

	double p1=sqr(q1),p2=sqr(q2),p3=sqr(q3),p4=sqr(q4);
	rotmx[0][0]=p1+p2-p3-p4;
	rotmx[1][1]=p1-p2+p3-p4;
	rotmx[2][2]=p1-p2-p3+p4;
	p1=q2*q3*2;p2=q1*q4*2;
	rotmx[0][1]=p1-p2; rotmx[1][0]=p1+p2;
	p1=q2*q4*2;p2=q1*q3*2;
	rotmx[0][2]=p1+p2; rotmx[2][0]=p1-p2;
	p1=q3*q4*2;p2=q1*q2*2;
	rotmx[1][2]=p1-p2; rotmx[2][1]=p1+p2;
}

// void MOLECULE::getRotMX()
// {
// 	double cf = std::cos(phi), sf = std::sin(phi), cq = std::cos(theta), sq = std::sin(theta), cy = std::cos(psi), sy = std::sin(psi);
// 	rotmx[0][0]=cf*cy - cq*sf*sy; rotmx[0][1]=-cq*cy*sf - cf*sy; rotmx[0][2]=  sf*sq;
// 	rotmx[1][0]=cy*sf + cf*cq*sy; rotmx[1][1]= cf*cq*cy - sf*sy; rotmx[1][2]= -cf*sq;
// 	rotmx[2][0]=sq*sy           ; rotmx[2][1]= cy*sq           ; rotmx[2][2]=  cq   ;
// }

// void PoincareBy(double **x_in,  const double &phi, const double &theta, const double &psi,
// double *dx, double **x_out, const int count)
// {
// 	double cf = std::cos(phi), sf = std::sin(phi), cq = std::cos(theta), sq = std::sin(theta), cy = std::cos(psi), sy = std::sin(psi);
// 	double rxy[3][3] = {{cf*cy - cq*sf*sy, -cq*cy*sf - cf*sy,  sf*sq}, 
// 	                    {cy*sf + cf*cq*sy,  cf*cq*cy - sf*sy, -cf*sq},
// 						{sq*sy           ,  cy*sq           ,  cq   }};
// 	if(dx)
// 	{
// 		for (int i = 0; i < count; ++i)
// 		{
//     	    x_out[i][0] = dx[0] + rxy[0][0]*x_in[i][0] + rxy[0][1]*x_in[i][1] + rxy[0][2]*x_in[i][2];
//     	    x_out[i][1] = dx[1] + rxy[1][0]*x_in[i][0] + rxy[1][1]*x_in[i][1] + rxy[1][2]*x_in[i][2];
//     	    x_out[i][2] = dx[2] + rxy[2][0]*x_in[i][0] + rxy[2][1]*x_in[i][1] + rxy[2][2]*x_in[i][2];
// 		}
// 	}
// 	else
// 	{
// 		for (int i = 0; i < count; ++i)
// 		{
//     	    x_out[i][0] = rxy[0][0]*x_in[i][0] + rxy[0][1]*x_in[i][1] + rxy[0][2]*x_in[i][2];
//     	    x_out[i][1] = rxy[1][0]*x_in[i][0] + rxy[1][1]*x_in[i][1] + rxy[1][2]*x_in[i][2];
//     	    x_out[i][2] = rxy[2][0]*x_in[i][0] + rxy[2][1]*x_in[i][1] + rxy[2][2]*x_in[i][2];
// 		}
// 	}
// }


// void AtomVelocity(double **x_in,  const double &phi, const double &theta, const double &psi,
// double *dx, double ***x_out, int count)
// {
// 	double cf = std::cos(phi), sf = std::sin(phi), cq = std::cos(theta), sq = std::sin(theta), cy = std::cos(psi), sy = std::sin(psi);
// 	double rxy[3][3][3] ={
//                        {{-cy*sf - cf*cq*sy,  -cf*cq*cy + sf*sy, cf*sq}, 
// 	                    { cf*cy - cq*sf*sy,  -cq*cy*sf - cf*sy, sf*sq},
// 						{ 0               ,   0               , 0    }},

//                        {{ sf*sq*sy,   cy*sf*sq,  cq*sf}, 
// 	                    {-cf*sq*sy,  -cf*cy*sq, -cf*cq},
// 						{ cq*sy   ,   cq*cy   , -sq   }},

//                       {{-cq*cy*sf - cf*sy,  -cf*cy + cq*sf*sy, 0}, 
// 	                   { cf*cq*cy - sf*sy,  -cy*sf - cf*cq*sy, 0},
// 					   { cy*sq           ,  -sq*sy           , 0}},
//     };

// 	for (int i = 0; i < count; ++i)
// 	{
//         for(int d=0;d<3;++d)
//         {
//             x_out[i][d][0] = dx[0] + rxy[d][0][0]*x_in[i][0] + rxy[d][0][1]*x_in[i][1] + rxy[d][0][2]*x_in[i][2];
//             x_out[i][d][1] = dx[1] + rxy[d][1][0]*x_in[i][0] + rxy[d][1][1]*x_in[i][1] + rxy[d][1][2]*x_in[i][2];
//             x_out[i][d][2] = dx[2] + rxy[d][2][0]*x_in[i][0] + rxy[d][2][1]*x_in[i][1] + rxy[d][2][2]*x_in[i][2];
//         }
// 	}
// }

void rotationMx(double *axis, const double &cc, double **mx)
{
	double ss=std::sqrt(1-cc*cc),ccc=1-cc;
	mx[0][0]=sqr(axis[0])*ccc+cc; mx[1][1]=sqr(axis[1])*ccc+cc; mx[2][2]=sqr(axis[2])*ccc+cc;
	mx[0][1]=mx[1][0]=axis[0]*axis[1]*ccc;
	mx[0][2]=mx[2][0]=axis[0]*axis[2]*ccc;
	mx[1][2]=mx[2][1]=axis[1]*axis[2]*ccc;
	double sx=ss*axis[0],sy=ss*axis[1],sz=ss*axis[2];
	mx[0][1]-=sz; mx[1][0]+=sz; 
	mx[0][2]+=sy; mx[2][0]-=sy; 
	mx[1][2]-=sx; mx[2][1]+=sx; 
}

void vecCross(double *v1,double *v2, double *vo)
{
	vo[0]=v1[1]*v2[2]-v1[2]*v2[1];
	vo[1]=v1[2]*v2[0]-v1[0]*v2[2];
	vo[2]=v1[0]*v2[1]-v1[1]*v2[0];
}
inline double vecDot(double *v1,double *v2)
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

inline double vecCosTheta(double *v1,double *v2)
{
	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/std::sqrt(vecDot(v1,v1)*vecDot(v2,v2));
}

inline void vecMul(double *vv, double sc)
{
	vv[0]*=sc;
	vv[1]*=sc;
	vv[2]*=sc;
}

void vecPerpComponent(double *vec, double *axis, double *v_out)
{
	// v - <v|a> a/|a|^2
	double fc=vecDot(vec,axis)/vecDot(axis,axis);
	v_out[0]=vec[0]-fc*axis[0];
	v_out[1]=vec[1]-fc*axis[1];
	v_out[2]=vec[2]-fc*axis[2];
}

// 1. generate rotation so that xi1 is rotated to xo1
// 2. apply rotation to xi2
// 3. generate rotation with xo1 as axis and rotated-xi2 is rotated to xo2
// 4. calculate matrix product of two rotation mx and convert to euler angles
void MOLECULE::computeEuler(double *xo1, double *xo2, double *xi1, double *xi2)
{
	double tht, axis[3], xit[3],xjt[3];

	// step 1
	vecCross(xi1,xo1,axis);
	tht=1/std::sqrt(vecDot(axis,axis)); vecMul(axis,tht);
	tht=vecCosTheta(xi1,xo1);
	rotationMx(axis,tht,tmpmx1);
	// step 2
	// lin_mxvec(tmpmx1[0],CblasNoTrans,xi2,3,3,0.0,xit);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,tmpmx1[0],3,xi2,1,0.0,xit,1);
	// step 3
	vecPerpComponent(xit,xo1,xit); // safe to call this with same in and out
	vecPerpComponent(xo2,xo1,xjt);
	vecCross(xit,xjt,axis);
	tht=1/std::sqrt(vecDot(axis,axis)); vecMul(axis,tht);	
	tht=vecCosTheta(xit,xjt);
	rotationMx(axis,tht,tmpmx2);
	// step 4
	// lin_mxmx(tmpmx2[0],CblasNoTrans,tmpmx1[0],CblasNoTrans,3,3,3,0.0,rotmx[0]);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.0,tmpmx2[0],3,tmpmx1[0],3,0.0,rotmx[0],3);
	theta=std::acos(rotmx[2][2]);
	psi=std::atan2(rotmx[2][0], rotmx[2][1]);
	phi=std::atan2(rotmx[0][2],-rotmx[1][2]);
}

inline void MOLECULE::writex()
{
	// PoincareBy(x_internal,phi,theta,psi,com,x+ibegin,nat);
	for(int ii=ibegin;ii<iend;++ii) memcpy(x[ii],com,3*sizeof(double));
	// lin_mxmx(x_internal[0],CblasNoTrans,rotmx[0],CblasTrans,nat,3,3,1.0,x[ibegin]);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,nat,3,3,1.0,x_internal[0],3,rotmx[0],3,1.0,x[ibegin],3);
}

void MOLECULE::rotateInternal()
{
	double temp;
	temp=innertia1; innertia1=innertia3; innertia3=temp;
	for(int i=0;i<nat;++i)
	{
		temp=x_internal[i][0]; x_internal[i][0]=-x_internal[i][2]; x_internal[i][2]=temp;
	}
}

void MOLECULE::init(int na, double **xi, int i_o, double *innert, int i1, int i2)
{
	nat=na;
	create2DArray(x_internal,na,3);
	for(int i=0;i<na;++i)
		for(int j=0;j<3;++j)
			x_internal[i][j]=xi[i][j];
	ibegin=i_o;
	iend=i_o+na;
	innertia1=innert[0]; innertia2=innert[1]; innertia3=innert[2];
	double xtemp;
	vcom[0]=vcom[1]=vcom[2]=angularl[0]=angularl[1]=angularl[2]=fcom[0]=fcom[1]=fcom[2]=torq[0]=torq[1]=torq[2]=0;
	r_ms=1.0/r_m[ibegin];
	com[0]=x[ibegin][0]*r_ms; com[1]=x[ibegin][1]*r_ms; com[2]=x[ibegin][2]*r_ms;
	for(int ii=ibegin+1;ii<iend;++ii)
	{
		r_ms+=1.0/r_m[ii];
		xtemp=x[ii][0]-x[ibegin][0]; xtemp-=std::nearbyint(xtemp*r_bl)*bl; com[0]+=((x[ii][0]=x[ibegin][0]+xtemp)/r_m[ii]);
		xtemp=x[ii][1]-x[ibegin][1]; xtemp-=std::nearbyint(xtemp*r_bl)*bl; com[1]+=((x[ii][1]=x[ibegin][1]+xtemp)/r_m[ii]);
		xtemp=x[ii][2]-x[ibegin][2]; xtemp-=std::nearbyint(xtemp*r_bl)*bl; com[2]+=((x[ii][2]=x[ibegin][2]+xtemp)/r_m[ii]);
	}
	r_ms=1.0/r_ms;
	com[0]*=r_ms; com[1]*=r_ms; com[2]*=r_ms;
	for(int ii=ibegin;ii<iend;++ii)
	{
		for(int jj=ibegin;jj<iend;++jj) special[ii][jj]=0;
		x[ii][0]-=com[0];
		x[ii][1]-=com[1];
		x[ii][2]-=com[2];
	}
	computeEuler(x[ibegin+i1],x[ibegin+i2],x_internal[i1],x_internal[i2]);
	// std::cerr<<"rotmx[2][2] = "<<rotmx[2][2]<<std::endl;
#ifdef MOL_AUTOROT
	if(rotmx[2][2]>ROT_THR || rotmx[2][2]<-ROT_THR) 
	{
		rotateInternal();
		computeEuler(x[ibegin+i1],x[ibegin+i2],x_internal[i1],x_internal[i2]);
		// std::cerr<<"x_internal rotated. new rotmx[2][2] = "<<rotmx[2][2]<<std::endl;
	}
#endif
	// std::cerr<<"Euler: "<<phi<<' '<<theta<<' '<<psi<<std::endl;
	com[0]-=std::nearbyint(com[0]*r_bl)*bl;
	com[1]-=std::nearbyint(com[1]*r_bl)*bl;
	com[2]-=std::nearbyint(com[2]*r_bl)*bl;
	writex();
	// std::cerr<<"det0 = "<<lin_det3(rotmx)<<std::endl;
}

// calculate COM force and torque
void MOLECULE::computeMolForce()
{
	double dx,dy,dz,fx,fy,fz;
	fcom[0]=fcom[1]=fcom[2]=torq[0]=torq[1]=torq[2]=0;
	for(int i=ibegin;i<iend;++i)
	{
		fx=f[i][0];fy=f[i][1];fz=f[i][2];
		fcom[0]+=fx;fcom[1]+=fy;fcom[2]+=fz;
		dx=x[i][0]-com[0]; dy=x[i][1]-com[1]; dz=x[i][2]-com[2];
		torq[0]+=(dy*fz-dz*fy);
		torq[1]+=(dz*fx-dx*fz);
		torq[2]+=(dx*fy-dy*fx);
	}

	fcom[0]*=(r_ms*accel); fcom[1]*=(r_ms*accel); fcom[2]*=(r_ms*accel);
	torq[0]*=accel; torq[1]*=accel; torq[2]*=accel;
}

void MOLECULE::updateX()
{
	double da;
	updateV();

	com[0]+=vcom[0]*dt;
	com[1]+=vcom[1]*dt;
	com[2]+=vcom[2]*dt;

	tmpmx4[0][0]=tmpmx4[1][1]=tmpmx4[2][2]=0;
	//init mx2 as rotmx
	memcpy(tmpmx2[0],rotmx[0],3*3*sizeof(double));

	for(int niter=0;niter<10;++niter)
	{
		//2 -> 1
		// printMatrix(tmpmx1,3,3);
		for(int i=0;i<3;++i)
		{
			tmpmx3[i][0]=rotmx[0][i]+tmpmx2[0][i];
			tmpmx3[i][1]=rotmx[1][i]+tmpmx2[1][i];
			tmpmx3[i][2]=rotmx[2][i]+tmpmx2[2][i];
		}
		// printMatrix(tmpmx3,3,3);
		// lin_mxvec(tmpmx3[0],CblasNoTrans,angularl,3,3,0.0,angularL); // angularL is used as omega here.
		cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,tmpmx3[0],3,angularl,1,0.0,angularL,1);
		angularL[0]/=innertia1; angularL[1]/=innertia2; angularL[2]/=innertia3; 
		tmpmx4[0][1]= angularL[2]; tmpmx4[0][2]=-angularL[1]; tmpmx4[1][2]= angularL[0]; 
		tmpmx4[1][0]=-angularL[2]; tmpmx4[2][0]= angularL[1]; tmpmx4[2][1]=-angularL[0]; 
		// printMatrix(tmpmx4,3,3);
		// lin_mxmx(tmpmx4[0],CblasNoTrans,tmpmx3[0],CblasNoTrans,3,3,3,0.0,tmpmx5[0]);
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.0,tmpmx4[0],3,tmpmx3[0],3,0.0,tmpmx5[0],3);
		// printMatrix(tmpmx1,3,3);
		for(int i=0;i<3;++i)
		{
			tmpmx1[i][0]=rotmx[i][0]+qdt*tmpmx5[0][i];
			tmpmx1[i][1]=rotmx[i][1]+qdt*tmpmx5[1][i];
			tmpmx1[i][2]=rotmx[i][2]+qdt*tmpmx5[2][i];
		}
		// printMatrix(tmpmx1,3,3,"mx1");

		//1 -> 2
		for(int i=0;i<3;++i)
		{
			tmpmx3[i][0]=rotmx[0][i]+tmpmx1[0][i];
			tmpmx3[i][1]=rotmx[1][i]+tmpmx1[1][i];
			tmpmx3[i][2]=rotmx[2][i]+tmpmx1[2][i];
		}
		// printMatrix(tmpmx3,3,3,"mx3");
		// lin_mxvec(tmpmx3[0],CblasNoTrans,angularl,3,3,0.0,angularL); // angularL is used as omega here.
		cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,tmpmx3[0],3,angularl,1,0.0,angularL,1);
		angularL[0]/=innertia1; angularL[1]/=innertia2; angularL[2]/=innertia3; 
		tmpmx4[0][1]= angularL[2]; tmpmx4[0][2]=-angularL[1]; tmpmx4[1][2]= angularL[0]; 
		tmpmx4[1][0]=-angularL[2]; tmpmx4[2][0]= angularL[1]; tmpmx4[2][1]=-angularL[0]; 
		// lin_mxmx(tmpmx4[0],CblasNoTrans,tmpmx3[0],CblasNoTrans,3,3,3,0.0,tmpmx5[0]);
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.0,tmpmx4[0],3,tmpmx3[0],3,0.0,tmpmx5[0],3);
		// printMatrix(rotmx,3,3,"mx2");
		for(int i=0;i<3;++i)
		{
			tmpmx2[i][0]=rotmx[i][0]+qdt*tmpmx5[0][i];
			tmpmx2[i][1]=rotmx[i][1]+qdt*tmpmx5[1][i];
			tmpmx2[i][2]=rotmx[i][2]+qdt*tmpmx5[2][i];
		}
		// printMatrix(tmpmx2,3,3,"mx2");

		// check convergence
		da=0;
		for(int i=0;i<3;++i)
		{
			da+=
			 sqr(tmpmx2[i][0]-tmpmx1[i][0])
			+sqr(tmpmx2[i][1]-tmpmx1[i][1])
			+sqr(tmpmx2[i][2]-tmpmx1[i][2]);
		}
		// std::cerr<<"err ~ "<<da<<std::endl;
		if(da<ITER_ACCU) break;
	}
	memcpy(rotmx[0],tmpmx2[0],3*3*sizeof(double));
	// the method here can ensure orthogonality itself, normalizing is actually optional
	nomalizeRotmx();
	writex();
}

void MOLECULE::updateV()
{
	vcom[0]+=fcom[0]*hdt;
	vcom[1]+=fcom[1]*hdt;
	vcom[2]+=fcom[2]*hdt;
	angularl[0]+=torq[0]*hdt;
	angularl[1]+=torq[1]*hdt;
	angularl[2]+=torq[2]*hdt;	
}

double MOLECULE::kineticErg()
{
	// lin_mxvec(rotmx[0],CblasTrans,angularl,3,3,0.0,angularL);
	cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.0,rotmx[0],3,angularl,1,0.0,angularL,1);
	return .5*((sqr(vcom[0])+sqr(vcom[1])+sqr(vcom[2]))/r_ms + sqr(angularL[0])/innertia1+sqr(angularL[1])/innertia2+sqr(angularL[2])/innertia3);
}

void MOLECULE::computeRotMK(double **tmx, double **mmx, double **smx)
{
	if(rotmx[2][2]>ROT_THR || rotmx[2][2]<-ROT_THR) std::cerr<<"MOLECULE::computeRotMMX -> Warning: potentially bad rotmx[2][2] "<<rotmx[2][2]<<std::endl;
	theta=std::acos(rotmx[2][2]);
	psi=std::atan2(rotmx[2][0], rotmx[2][1]);
	phi=std::atan2(rotmx[0][2],-rotmx[1][2]);

	double cf = std::cos(phi), sf = std::sin(phi), cq = std::cos(theta), sq = std::sin(theta), cy = std::cos(psi), sy = std::sin(psi);

	double cfcq=cf*cq, cfcy=cf*cy, cqcy=cq*cy;
	double sfsq=sf*sq, sfsy=sf*sy, sqsy=sq*sy;
	double cfsq=cf*sq, cfsy=cf*sy, cqsy=cq*sy;
	double cqsf=cq*sf, cysf=cy*sf, cysq=cy*sq;
	double cfcqcy=cf*cq*cy;
	double sfsqsy=sf*sq*sy;
	double cfcqsy=cf*cq*sy, cfcysq=cf*cy*sq, cqcysf=cq*cy*sf;
	double cfsqsy=cf*sq*sy, cqsfsy=cq*sf*sy, cysfsq=cy*sf*sq;

	mmx[0][0]=innertia3*sqr(cq)+sqr(sq)*(innertia2*sqr(cy)+innertia1*sqr(sy));
	mmx[1][1]=.5*(innertia1+innertia2+(innertia1-innertia2)*(sqr(cy)-sqr(sy)));
	mmx[2][2]=innertia3;
	mmx[0][1]=(innertia1-innertia2)*sqsy*cy;
	mmx[0][2]=innertia3*cq;
	mmx[1][2]=mmx[2][1]=mmx[1][0]=mmx[2][0]=0; // also set lower triangle to zero

	int tnat=3*nat;
	double *tmmx=tmx[0];
	for(int i=2;i<tnat;i+=3) tmmx[i]=0; //zero tmx[0][3i+2] for simplification of calculation

	// phi phi
    tmpmx3[0][0]= -cfcy + cqsfsy; tmpmx3[1][0]=  cqcysf + cfsy; tmpmx3[2][0]= -sfsq;
	tmpmx3[0][1]= -cysf - cfcqsy; tmpmx3[1][1]= -cfcqcy + sfsy; tmpmx3[2][1]=  cfsq;
	// tmpmx3[0][2]=  0            ; tmpmx3[1][2]=   0           ; tmpmx3[2][2]= 0    ;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,2,3,1.0,x_internal[0],3,tmpmx3[0],3,0.0,tmmx,3);
	smx[0][0]=cblas_ddot(tnat,tmmx,1,f[ibegin],1);

	// psi psi
    //tmpmx3[0][0]= -cfcy + cqsfsy; tmpmx3[1][0]=  cqcysf + cfsy; tmpmx3[2][0]= 0;
	//tmpmx3[0][1]= -cysf - cfcqsy; tmpmx3[1][1]= -cfcqcy + sfsy; tmpmx3[2][1]= 0;
	tmpmx3[0][2]=  -sqsy        ; tmpmx3[1][2]=  -cysq        ;// tmpmx3[2][2]= 0;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,3,2,1.0,x_internal[0],3,tmpmx3[0],3,0.0,tmx[1],3);
	smx[2][2]=cblas_ddot(tnat,tmx[1],1,f[ibegin],1);

	// phi psi
    tmpmx5[0][0]= -cfcqcy + sfsy; tmpmx5[1][0]=  cysf + cfcqsy; //tmpmx5[2][0]= 0;
	tmpmx5[0][1]= -cqcysf - cfsy; tmpmx5[1][1]= -cfcy + cqsfsy; //tmpmx5[2][1]= 0;
	// tmpmx5[0][2]=  0            ; tmpmx5[1][2]=  0            ; tmpmx5[2][2]= 0;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,2,2,1.0,x_internal[0],3,tmpmx5[0],3,0.0,tmmx,3);
	smx[0][2]=smx[2][0]=cblas_ddot(tnat,tmmx,1,f[ibegin],1);

	// theta theta
    tmpmx3[0][0]=   cqsfsy; tmpmx3[1][0]=  cqcysf; tmpmx3[2][0]= -sfsq;
	tmpmx3[0][1]=  -cfcqsy; tmpmx3[1][1]= -cfcqcy; tmpmx3[2][1]= cfsq;
	tmpmx3[0][2]=  -sqsy  ; tmpmx3[1][2]=  -cysq ; tmpmx3[2][2]= -cq;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,3,3,1.0,x_internal[0],3,tmpmx3[0],3,0.0,tmx[1],3);
	smx[1][1]=cblas_ddot(tnat,tmx[1],1,f[ibegin],1);

	// phi theta
    tmpmx4[0][0]= cfsqsy; tmpmx4[1][0]= cfcysq; tmpmx4[2][0]= cfcq;
	tmpmx4[0][1]= sfsqsy; tmpmx4[1][1]= cysfsq; tmpmx4[2][1]= cqsf;
	// tmpmx4[0][2]=  0    ; tmpmx4[1][2]=   0   ; tmpmx4[2][2]= 0   ;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,2,3,1.0,x_internal[0],3,tmpmx4[0],3,0.0,tmmx,3);
	smx[0][1]=smx[1][0]=cblas_ddot(tnat,tmmx,1,f[ibegin],1);

	// theta psi
    tmpmx5[0][0]=  cysfsq; tmpmx5[1][0]= -sfsqsy;// tmpmx5[2][0]= 0;
	tmpmx5[0][1]= -cfcysq; tmpmx5[1][1]=  cfsqsy;// tmpmx5[2][1]= 0;
	tmpmx5[0][2]=  cqcy  ; tmpmx5[1][2]= -cqsy  ;// tmpmx5[2][2]= 0;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,3,2,1.0,x_internal[0],3,tmpmx5[0],3,0.0,tmx[1],3);
	smx[1][2]=smx[2][1]=cblas_ddot(tnat,tmx[1],1,f[ibegin],1);

	// phi
    tmpmx3[0][0]= -cysf - cfcqsy; tmpmx3[1][0]=  -cfcqcy + sfsy; tmpmx3[2][0]= cfsq;
	tmpmx3[0][1]=  cfcy - cqsfsy; tmpmx3[1][1]=  -cqcysf - cfsy; tmpmx3[2][1]= sfsq;
	// tmpmx3[0][2]=  0            ; tmpmx3[1][2]=   0            ; tmpmx3[2][2]= 0   ;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,2,3,1.0,x_internal[0],3,tmpmx3[0],3,0.0,tmmx,3);

	// theta
	tmpmx4[0][0]=  sfsqsy        ; tmpmx4[1][0]=  cysfsq         ; tmpmx4[2][0]=  cqsf; 
	tmpmx4[0][1]= -cfsqsy        ; tmpmx4[1][1]= -cfcysq         ; tmpmx4[2][1]= -cfcq;
	tmpmx4[0][2]=  cqsy          ; tmpmx4[1][2]=  cqcy           ; tmpmx4[2][2]= -sq  ;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,3,3,1.0,x_internal[0],3,tmpmx4[0],3,0.0,tmx[1],3);

	// psi
    tmpmx5[0][0]= -cqcysf - cfsy; tmpmx5[1][0]=  -cfcy + cqsfsy;// tmpmx5[2][0]= 0;
	tmpmx5[0][1]=  cfcqcy - sfsy; tmpmx5[1][1]=  -cysf - cfcqsy;// tmpmx5[2][1]= 0;
	tmpmx5[0][2]=  cysq         ; tmpmx5[1][2]=  -sqsy         ;// tmpmx5[2][2]= 0;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nat,3,2,1.0,x_internal[0],3,tmpmx5[0],3,0.0,tmx[2],3);

}


// void MOLECULE::computeMolForce()
// {
// 	double dx,dy,dz,fx,fy,fz,tx,ty,tz;
// 	fcom[0]=fcom[1]=fcom[2]=tx=ty=tz=0;
// 	for(int i=ibegin;i<iend;++i)
// 	{
// 		fx=f[i][0];fy=f[i][1];fz=f[i][2];
// 		fcom[0]+=fx;fcom[1]+=fy;fcom[2]+=fz;
// 		dx=x[i][0]-com[0]; dy=x[i][1]-com[1]; dz=x[i][2]-com[2];
// 		tx+=(dy*fz-dz*fy);
// 		ty+=(dz*fx-dx*fz);
// 		tz+=(dx*fy-dy*fx);
// 	}
// 	getAxis();
// 	torq1=tx*rotmx[0][0]+ty*rotmx[1][0]+tz*rotmx[2][0];
// 	torq2=tx*rotmx[0][1]+ty*rotmx[1][1]+tz*rotmx[2][1];
// 	torq3=tx*rotmx[0][2]+ty*rotmx[1][2]+tz*rotmx[2][2];

// 	fcom[0]*=(r_ms*accel); fcom[1]*=(r_ms*accel); fcom[2]*=(r_ms*accel);
// 	dw1=(torq1+(innertia2-innertia3)*omega2*omega3)/innertia1*accel;
// 	dw2=(torq2+(innertia3-innertia1)*omega3*omega1)/innertia2*accel;
// 	dw3=(torq3+(innertia1-innertia2)*omega1*omega2)/innertia3*accel;
// }

// void MOLECULE::computeV_euler()
// {
// 	double cy=std::cos(psi),sy=std::sin(psi);
// 	v_phi=(omega1*sy+omega2*cy)/std::sin(theta);
// 	v_theta=(omega1*cy-omega2*sy);
// 	v_psi=omega3-v_phi*std::cos(theta);
// }

// void MOLECULE::updateX()
// {
// 	updateV();

// 	com[0]+=vcom[0]*dt;
// 	com[1]+=vcom[1]*dt;
// 	com[2]+=vcom[2]*dt;
// 	computeV_euler();
// 	phi+=(v1_phi=v_phi)*dt;
// 	theta+=(v1_theta=v_theta)*dt;
// 	psi+=(v1_psi=v_psi)*dt;
// 	writex();
// }

// void MOLECULE::updateV()
// {
// 	vcom[0]+=fcom[0]*hdt;
// 	vcom[1]+=fcom[1]*hdt;
// 	vcom[2]+=fcom[2]*hdt;
// 	omega1+=dw1*hdt;
// 	omega2+=dw2*hdt;
// 	omega3+=dw3*hdt;	
// }

// double MOLECULE::kineticErg()
// {
// 	return .5*((sqr(vcom[0])+sqr(vcom[1])+sqr(vcom[2]))/r_ms + innertia1*sqr(omega1)+innertia2*sqr(omega2)+innertia3*sqr(omega3));
// }