#include "src/includes.hpp"
#include <cmath>



double U(double x, double y)
{
    return 0.0;
}

double V(double x, double y)
{
    return 0.0;
}

double eta(double x, double y)
{
    return (exp(-8*((x)*(x)+y*y)));
}

double Depth(double x, double y)
{
    return 0.0;
}
int main()
{
    unsigned    Nex =   20;
    unsigned    Ney =   20;
    unsigned    N   =   4;
    double L_start  =   -1;
    double L_end    =   1;
    double H_start  =   -1;
    double H_end    =   1;
    ShallowWater q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setDepth(Depth);
    q.setInitialConditions(eta,U,V);
    q.setSolver(1e-4,400);
    q.solve();
    q.plotSolution("t=0.001s");

    return 0;
}
