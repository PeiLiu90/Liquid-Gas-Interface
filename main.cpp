#include "Liquid-Liquid.h"
#include "parameter.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "time.h"
#include <stdlib.h>
using namespace std;

int main()
{
    vector2d<double,N1,N1> sig;
    vector2d<double,N1,N1> eps;
    vector2d<double,N1,N1> alp;
    Load(sig,eps,alp);
    //Liquid_liquid<N1,M1> LL(sig,eps,alp, l,rho);
    RHNC<N1,M1-1> LJ(sig,eps,alp,l,rho);
    LJ.ComputeCorrelationFunction();

    return 0;
}


