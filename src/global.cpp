
#include "global.h"
#include "memory.h"
#include "molecule.h"
#include <cstring>
double **x, **v, **f; // position, velocity, force
double *q, *r_m; // charge, inv_mass
double **special; //special bond
int *typ;
int natom=0, ntype;//, ntotal, x_capacity;
// double **g_dx;
// int *g_orig;
int nmol=0;
MOLECULE *molecule;
// we assume that lo=-h_bl and hi=h_bl;
double bl=0., r_bl, h_bl; // box_length, inverse_bl, half_bl

const double pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286;
const double two_pi=pi*2, four_pi=pi*4;

const double electricK=332.06371, epsilon_0=1/(four_pi*electricK);
const double accel=4.184e-4;

void initialize(int ntp, const char *xyzfile, const char *elms[], const double *r_mass, const double *charges)
{
    ntype=ntp;
    char buf[256];
    char elm[256];
    int iel;
    FILE *fp=fopen(xyzfile,"r");

    fgets(buf,256,fp);
    sscanf(buf,"%d",&natom);
    std::cerr<<"natom = "<<natom<<std::endl;
    if (natom==0) throw std::runtime_error("initialize: error reading natom");

    fgets(buf,256,fp);
    if(sscanf(buf,"%*s %*s %*s %*s %lf %lf",&h_bl,&r_bl) != 2) throw std::runtime_error("initialize: error reading box length");
    set_bl(r_bl-h_bl);
    std::cerr<<"bl = "<<bl<<std::endl;

    create1DArray(typ,natom);
    create1DArray(q,natom);
    create1DArray(r_m,natom);
    create2DArray(special,natom,natom);
    create2DArray(x,natom,3);
    create2DArray(v,natom,3);
    create2DArray(f,natom,3);
    for(int idd=0;idd<natom;++idd)
    {
        for(int gx=0;gx<natom;++gx) special[idd][gx]=1;
        fgets(buf,256,fp);
        if(sscanf(buf,"%s %lf %lf %lf",elm,x[idd],x[idd]+1,x[idd]+2) != 4) throw std::runtime_error("initialize: error reading atoms");
        for(iel=0;iel<ntype;++iel)
            if(strcmp(elm,elms[iel]) == 0) break;
        if(iel==ntype) throw std::runtime_error("initialize: unknown element");
        typ[idd]=iel;
        r_m[idd]=r_mass[iel];
        q[idd]=charges[iel];
    }
    
    fclose(fp);
}

int readMolecule(const char*xyzfile, const char *elms[], const double *r_mass, double **&xin, double *innertia)
{
    char buf[256];
    char elm[256];
    int na,iel;
    FILE *fp=fopen(xyzfile,"r");

    fgets(buf,256,fp);
    sscanf(buf,"%d",&na);
    if (na==0) throw std::runtime_error("readMolecule: error reading natom");

    fgets(buf,256,fp);
    // if(sscanf(buf,"%*s %lf %lf %lf",innertia,innertia+1,innertia+2) != 2) throw std::runtime_error("readMolecule: error reading innertia");
    // std::cerr<<"innertia = "<<innertia[0]<<' '<<innertia[1]<<' '<<innertia[2]<<std::endl;

    create2DArray(xin,na,3);
    double *ms=new double[na];
    double com[3]={0,0,0},totm=0;
    for(int idd=0;idd<na;++idd)
    {
        fgets(buf,256,fp);
        if(sscanf(buf,"%s %lf %lf %lf",elm, xin[idd],xin[idd]+1,xin[idd]+2) != 4) throw std::runtime_error("readMolecule: error reading atoms");
        for(iel=0;iel<ntype;++iel)
            if(strcmp(elm,elms[iel]) == 0) break;
        if(iel==ntype) throw std::runtime_error("initialize: unknown element");
        totm+=(ms[idd]=1./r_mass[iel]);
        com[0]+=(xin[idd][0]*ms[idd]);
        com[1]+=(xin[idd][1]*ms[idd]);
        com[2]+=(xin[idd][2]*ms[idd]);
    }
    fclose(fp);
    com[0]/=totm; com[1]/=totm; com[2]/=totm;
    innertia[0]=innertia[1]=innertia[2]=0;
    double x2[3];
    for(int idd=0;idd<na;++idd)
    {
        x2[0]=sqr(xin[idd][0]-=com[0]);
        x2[1]=sqr(xin[idd][1]-=com[1]);
        x2[2]=sqr(xin[idd][2]-=com[2]);
        innertia[0]+=ms[idd]*(x2[1]+x2[2]);
        innertia[1]+=ms[idd]*(x2[2]+x2[0]);
        innertia[2]+=ms[idd]*(x2[0]+x2[1]);
        // std::cerr<<xin[idd][0]<<' '<<xin[idd][1]<<' '<<xin[idd][2]<<'\n';
    }
    delete[] ms;
    std::cerr<<"readMolecule: N = "<<na<<"; M = "<< totm<<"; I[] = "<<innertia[0]<<' '<<innertia[1]<<' '<<innertia[2]<<std::endl;
    return na;
}

//i1 and i2 is used to calculate the orientation of initial condition, so x[i1], x[i2] and COM must not be colinear.
void groupMol(int na, int nm, int ibegin, double **xi, double *innertia, int i1,int i2,double dt)
{
    nmol+=nm;
    molecule=new MOLECULE[nm];
    for(int i=0;i<nm;++i)
    {
        molecule[i].init(na,xi,ibegin,innertia,i1,i2);
        molecule[i].setDt(dt);
        ibegin+=na;
    }
    std::cerr<<"created "<<nm<<" molecules\n";
}

void printMatrix(double **mx,int m,int n,const char* info, char fmt)
{
	std::cerr<<info<<std::endl;
    switch (fmt)
    {
    case 'm':
        std::cerr<<'{'<<std::endl;
    	for(int i=0;i<m;++i)
	    {
            std::cerr<<'{';
	    	for(int j=0;j<n;++j)
            {
	    		std::cerr<<mx[i][j];
                if (j+1<n) std::cerr<<',';
            }
            if(i+1<m) std::cerr<<"},\n";
            else std::cerr<<'}'<<std::endl;
	    }
        std::cerr<<'}'<<std::endl;
        break;
    case 'p':
        std::cerr<<'['<<std::endl;
    	for(int i=0;i<m;++i)
	    {
            std::cerr<<'[';
	    	for(int j=0;j<n;++j)
            {
	    		std::cerr<<mx[i][j];
                if (j+1<n) std::cerr<<',';
            }
            if(i+1<m) std::cerr<<"],\n";
            else std::cerr<<']'<<std::endl;
	    }
        std::cerr<<']'<<std::endl;
        break;

        break;
    default:
    	for(int i=0;i<m;++i)
	    {
	    	for(int j=0;j<n;++j)
	    		std::cerr<<mx[i][j]<<' ';
	    	std::cerr<<std::endl;
	    }
        break;
    }
	std::cerr<<std::endl;
}

void writeArrayBin(double *arr,int n,const char* fname)
{
    FILE *fp=fopen(fname,"wb");
    fwrite(arr,sizeof(double),n,fp);
    fclose(fp);
}
