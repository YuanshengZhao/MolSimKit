#pragma once

class NEIGHLIST
{
protected:
    double cutsq;
    double **x_prev;
    double dx2;
    int capacity;
    // void refold();
public:
    int *num_neigh;
    int **nei_list;
    double ***nei_dx;
    int nbuild;
    NEIGHLIST(double cutoff, double skin, int cap);
    ~NEIGHLIST();
    template <bool full> void build(bool check=true); // check displacement before building do not build if not needed
    bool isvalid();
};

