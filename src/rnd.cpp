#include "rnd.h"
#include <random>
#include <iostream> 

std::mt19937 rnd_mt;
const auto rnd_max=rnd_mt.max();

void randomSetSeed(unsigned seed)
{
    if(seed>0) rnd_mt.seed(seed);
    else
    {
        std::random_device rd;
        seed=rd();
        std::cerr<<"randomSetSeed: "<<seed<<std::endl;
        rnd_mt.seed(seed);
    }
}

unsigned randIntStrict(unsigned n)
{
    auto minx=rnd_max%n, rd=rnd_mt();
    while(rd<=minx) rd=rnd_mt();
    return rd%n;
}
