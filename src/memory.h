#pragma once
#include<iostream>

template <typename numtyp>
numtyp* create1DArray(numtyp *&arr, int i);

template <typename numtyp>
void destroy1DArray(numtyp *&arr);

template <typename numtyp>
numtyp** create2DArray(numtyp **&arr, int i, int j);

template <typename numtyp>
void destroy2DArray(numtyp **&arr);

template <typename numtyp>
numtyp*** create3DArray(numtyp ***&arr, int i, int j, int k);

template <typename numtyp>
void destroy3DArray(numtyp ***&arr);