#include <iostream>
#include "memory.h"

template double*   create1DArray<double>(double *&,   int);
template double**  create2DArray<double>(double **&,  int, int);
template double*** create3DArray<double>(double ***&, int, int, int);

template int*   create1DArray<int>(int *&,   int);
template int**  create2DArray<int>(int **&,  int, int);
template int*** create3DArray<int>(int ***&, int, int, int);

template void destroy1DArray<double>(double *&);
template void destroy2DArray<double>(double **&);
template void destroy3DArray<double>(double ***&);

template void destroy1DArray<int>(int *&);
template void destroy2DArray<int>(int **&);
template void destroy3DArray<int>(int ***&);

template <typename numtyp>
numtyp* create1DArray(numtyp *&arr, int i)
{
    arr=new numtyp[i];
    return arr;
}

template <typename numtyp>
void destroy1DArray(numtyp *&arr)
{
    delete[] arr;
    arr=nullptr;
}

template <typename numtyp>
numtyp** create2DArray(numtyp **&arr, int i, int j)
{
    numtyp *dat=new numtyp[i*j];
    arr=new numtyp*[i];
    arr[0]=dat;
    for(int _i=1;_i<i;++_i) arr[_i]=arr[_i-1]+j;
    return arr;
}

template <typename numtyp>
void destroy2DArray(numtyp **&arr)
{
    delete[] arr[0];
    delete[] arr;
    arr=nullptr;
}


template <typename numtyp>
numtyp*** create3DArray(numtyp ***&arr, int i, int j, int k)
{
    numtyp *dat=new numtyp[i*j*k];
    numtyp **ptr=new numtyp*[i*j];
    arr=new numtyp**[i];
    ptr[0]=dat;
    int ij=i*j;
    for(int _i=1;_i<ij;++_i) ptr[_i]=ptr[_i-1]+k;
    arr[0]=ptr;
    for(int _i=1;_i<i;++_i) arr[_i]=arr[_i-1]+j;
    return arr;
}

template <typename numtyp>
void destroy3DArray(numtyp ***&arr)
{
    delete[] arr[0][0];
    delete[] arr[0];
    delete[] arr;
    arr=nullptr;
}