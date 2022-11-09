#pragma once

class PHONON
{
public:
    double **KMX; // K matrix on atom cord
    double **KMX_mol; // K matrix on molecule cord
    double ***MMX; // rotational M matrix for each molecule
    double *rr_mass; // inverse mass of molecules
    double ***transferMX; // jacobian D(atom cord)/D(euler angle)
    double *w2; // eigenvalues
    double **dr; // displacement vec in each mode
    PHONON();
    ~PHONON();
    // must build KMX outside first!, currently only support 1 type of molelcule
    // also must compute force outside first!
    void buildKMX_mol(); 
    void eigenKMX_mol();
    void symKMX();
};