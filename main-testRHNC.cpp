#include "RHNC.h"
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

    Load(sig,eps);

    RHNC<N1,M1> LJ(eps,sig,l,rho);

    LJ.ComputeCorrelationFunction();

    ofstream fout("DCF.txt");
    ofstream gout("RDF2.txt");
    ofstream bout("bridge.txt");
    for(int i=0;i<M1;i++)
    {
        fout<<LJ._r[i]<<"  ";
        gout<<LJ._r[i]<<"  ";
        bout<<LJ._r[i]<<"  ";
        for(int j=0;j<N1;j++)
        {
            for(int k=0;k<N1;k++)
            {
                fout<<LJ._C(j,k,i)<<"  ";
                gout<<LJ._H(j,k,i)+1.<<"  ";
                bout<<LJ.HardSphere._B(j,k,i)<<"  ";
            }
        }
        fout<<endl;
        gout<<endl;
        bout<<endl;
    }
    fout.close();
    gout.close();
    bout.close();
    return 0;
}


