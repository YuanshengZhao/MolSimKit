#pragma once
#include "neighlist.h"

class VERLET
{
private:
    double dt,hdt;
    double **f_prev;
public:
    VERLET(double ddt);
    ~VERLET();
    //usage: updateX => calc force => convert f->a => updateV
    void setDt(double ddt);
    void updateX();
    void updateV();
    void f2a();
    void scaleV(double sc);
};
