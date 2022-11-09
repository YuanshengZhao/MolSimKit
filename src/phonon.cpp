#include "phonon.h"
#include "global.h"
#include "memory.h"
#include "molecule.h"
#include "linalg.h"
#include "../lib/OpenBLAS/include/lapacke.h"

PHONON::PHONON():
KMX(nullptr),KMX_mol(nullptr),MMX(nullptr),rr_mass(nullptr),transferMX(nullptr),w2(nullptr),dr(nullptr)
{
    create2DArray(KMX,3*natom,3*natom);
}

PHONON::~PHONON()
{
    destroy2DArray(KMX);
    if(KMX_mol) destroy2DArray(KMX_mol);
    if(MMX) destroy3DArray(MMX);
    if(rr_mass) destroy1DArray(rr_mass);
    if(transferMX) destroy3DArray(transferMX);
    if(w2) destroy1DArray(w2);
    if(dr) destroy2DArray(dr);
}

void PHONON::symKMX()
{
    int tam=3*natom;
    for(int i=0;i<tam;++i)
        for(int j=i+1;j<tam;++j)
            KMX[j][i]=KMX[i][j];
}

void PHONON::buildKMX_mol()
{
    int s_i,tnai,tnai_,tmpi1,tmpi2;
    const int na=natom/nmol;
    const int tatm=3*natom,s_nm=6*nmol,tna=3*na;
    if (!KMX_mol) create2DArray(KMX_mol,s_nm,s_nm);
    if (!MMX)     create3DArray(MMX,nmol,3,3);
    if (!transferMX) create3DArray(transferMX,nmol,3,tna);
    double **tempmx1,*tpv1,*tpv2, ***kmx2;
    create3DArray(kmx2,nmol,3,3);
    create2DArray(tempmx1,s_nm,tatm);
    for(int i=0;i<nmol;++i)
    {
        if(na!=molecule[i].nat) throw std::runtime_error("PHONON::buildKMX_mol -> wrong n_atoms per mol");
        molecule[i].computeRotMK(transferMX[i],MMX[i],kmx2[i]);
    }

    //zero tempmx1 KMX_mol
    for(int i=0;i<s_nm;++i)
    {
        tpv1=tempmx1[i];
        for(int j=0;j<tatm;++j)
            tpv1[j]=0;
        tpv1=KMX_mol[i];
        for(int j=0;j<s_nm;++j)
            tpv1[j]=0;
    }

    /*
    The matrix to diagonolize is U . tr . K . tr^T . U^T  + U . kmx2 . U^T
    where U is the cholesky of M^{-1}, 
    tr is the transfer matrix from atom cord to molecule cord
    clearly it is more convenient to compute U . tr first
    */

    // step1.1: invert MMX and Cholesky
    // step1.2: compute rotational part of U . tr (transiational part is trivial) 
    double tempms;
    if (!rr_mass) create1DArray(rr_mass,nmol);
    for(int i=0;i<nmol;++i)
    {
        rr_mass[i]=std::sqrt(molecule[i].r_ms);
        // printMatrix(MMX[i],3,3,"M");
        if(lin_inverse_sym_pos(MMX[i],3)) throw std::runtime_error("PHONON::buildKMX_mol -> inverse M failed");
        // printMatrix(MMX[i],3,3,"M");
        if(LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',3,MMX[i][0],3)) throw std::runtime_error("PHONON::buildKMX_mol -> Cholesky M^{-1} failed");
        // printMatrix(MMX[i],3,3,"U");
        // printMatrix(transferMX[i],3,tna,"tr");
        cblas_dtrmm(CblasRowMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,3,tna,1.0,MMX[i][0],3,transferMX[i][0],tna); // now tr is overwritten by U . tr
        // printMatrix(transferMX[i],3,tna,"U . tr");
    }

    // step2: compute tempmx1 = transfer . KMX
    for(int i=0;i<nmol;++i)
    {
        /* 
        translational: accumulate KMX to tempmx
        for _i = 0 : n
        tempmx1[6i:6i+3][0:3N] += KMX[3ni+3_i:3ni+3_i+3][0:3N] * rr_mass[i]
        */
        s_i=6*i; tnai=tna*i;
        for(int _i=0;_i<tna;_i+=3)
        {
            tnai_=tnai+_i;
            for(int _j=0;_j<3;++_j)
            {
                tpv1=tempmx1[s_i+_j];
                tpv2=KMX[tnai_+_j];
                for(int _k=0;_k<tatm;++_k)
                    tpv1[_k]+=tpv2[_k];
            }
        }
        tempms=rr_mass[i];
        for(int _j=0;_j<3;++_j)
        {
            tpv1=tempmx1[s_i+_j];
            for(int _k=0;_k<tatm;++_k)
                tpv1[_k]*=tempms;
        }
        //rotational: matrix product
        //tempmx1[6i+3:6i+6][0:3N] = transfer[i] . KMX[3ni:3ni+3n][0:3N]
        // lin_mxmx(transferMX[i][0],CblasNoTrans,KMX[tnai],CblasNoTrans,3,tatm,tna,0.0,tempmx1[s_i+3]);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,tatm,tna,1.0,transferMX[i][0],tna,KMX[tnai],tatm,0.0,tempmx1[s_i+3],tatm);
    }

    // step3: compute KMX_mol = tempmx1 . transfer^T
    for(int i=0;i<nmol;++i)
    {
        /* 
        translational: accumulate tempmx to KMX_mol
        for _i = 0 : n
        KMX_mol[0:6M][6i:6i+3] += tempmx[0:6M][3ni+3_i:3ni+3_i+3]
        */
        s_i=6*i; tnai=tna*i;
        for(int _i=0;_i<tna;_i+=3)
        {
            tnai_=tnai+_i;
            for(int _k=0;_k<3;++_k)
            {
                tmpi1=_k+s_i;
                tmpi2=_k+tnai_;
                for(int _j=0;_j<s_nm;++_j)
                    KMX_mol[_j][tmpi1]+=tempmx1[_j][tmpi2];
            }
        }
        tempms=rr_mass[i];
        for(int _k=0;_k<3;++_k)
        {
            tmpi1=_k+s_i;
            for(int _j=0;_j<s_nm;++_j)
                KMX_mol[_j][tmpi1]*=tempms;
        }

        //rotational: matrix product
        //KMX[0:6M][6i+3:6i+6] = tempmx[0:6M][3ni:3ni+3n] . transfer[i]^T
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,s_nm,3,tna,1.0,tempmx1[0]+tnai,tatm,transferMX[i][0],tna,0.0,KMX_mol[0]+s_i+3,s_nm);
    }

    // step4: compute U . kmx2 . U^T and subtract from KMX_mol (note as force is used, there is a minus sign here)
    for(int i=0;i<nmol;++i)
    {
        cblas_dtrmm(CblasRowMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,3,3,1.0,MMX[i][0],3,kmx2[i][0],3);
        cblas_dtrmm(CblasRowMajor,CblasRight,CblasUpper,CblasTrans,CblasNonUnit,3,3,1.0,MMX[i][0],3,kmx2[i][0],3);
        // KMX_mol[6i+3:6i+6][6i+3:6i+6] -= kmx2[i]
        s_i=6*i+3;tnai=s_i+1;tnai_=tnai+1;
        KMX_mol[s_i][s_i]-=kmx2[i][0][0];
        KMX_mol[tnai][tnai]-=kmx2[i][1][1];
        KMX_mol[tnai_][tnai_]-=kmx2[i][2][2];
        KMX_mol[s_i][tnai]=(KMX_mol[tnai][s_i]-=kmx2[i][0][1]);
        KMX_mol[s_i][tnai_]=(KMX_mol[tnai_][s_i]-=kmx2[i][0][2]);
        KMX_mol[tnai][tnai_]=(KMX_mol[tnai_][tnai]-=kmx2[i][1][2]);
    }


    // printMatrix(KMX_mol,s_nm,s_nm,"KMol",'m');

    destroy2DArray(tempmx1);
    destroy3DArray(kmx2);
}

void PHONON::eigenKMX_mol()
{
    int s_i,tnai,tnai_,tpi1,tpi2;
    const int na=natom/nmol;
    const int tatm=3*natom,s_nm=6*nmol,tna=3*na;
    if (!w2) create1DArray(w2,s_nm);
    if (!dr) create2DArray(dr,s_nm,tatm);

    // step 1: compute eigen val and vec
    // now KMX_mol contains vecs as row vectors
    if (lin_eigensystem_sym(KMX_mol,s_nm,w2)) throw std::runtime_error("PHONON::eigenKMX_mol -> cannot compute eigen vals and vecs");

    // step 2: compute U . transfer
    // this is already done in building KMX_mol and stored in transferMX

    // step 3: compute J(now = KMX_mol) . transfer
    double rms;
    double *tpv1;
    for(int i=0;i<nmol;++i)
    {
        // translational:
        // for j = 0 to na
        // dr[:][3 na i + 3 j: 3 na i + 3 j + 3] = rr_mass[i] * (KMX_mol[6 i : 6 i + 3][:])^T
        tnai=tna*i;
        s_i=6*i;
        rms=rr_mass[i];
        // first atom is computed directly
        for(int _j=0;_j<3;++_j)
        {
            tpv1=KMX_mol[s_i+_j];
            tnai_=tnai+_j;
            for(int _i=0;_i<s_nm;++_i)
            {
                dr[_i][tnai_]=rms*tpv1[_i];
            }
        }
        // others are copied from first atom
        for(int j=3;j<tna;j+=3)
        {
            tnai_=tnai+j;
            for(int _j=0;_j<3;++_j)
            {
                tpi1=tnai_+_j;
                tpi2=tnai+_j;
                for(int _i=0;_i<s_nm;++_i)
                {
                    dr[_i][tpi1]=dr[_i][tpi2];
                }
            }
        }

        // rotational: 
        // dr[:][3 na i: 3 na i + 3 na] += (KMX_mol[6 i + 3: 6 i + 6][:])^T . transferMX[i]
        // Careful! it is += here!
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,s_nm,tna,3,1.0,KMX_mol[s_i+3],s_nm,transferMX[i][0],tna,1.0,dr[0]+tnai,tatm);
    }

}