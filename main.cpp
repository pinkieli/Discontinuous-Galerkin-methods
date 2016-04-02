#include "src/includes.hpp"
#include <cmath>

#define U0  1.0

double U(double x, double y)
{
    return y;
}

double V(double x, double y)
{
    return (-x);
}

double initialConditions(double x, double y)
{
    return (exp(-32*((x+0.5)*(x+0.5)+y*y)));
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
    Field q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setVelocity(U,V);
    q.setInitialConditions(initialConditions);
    q.setSolver(1e-3,1570);
    q.solve();
    q.plotSolution("t=1.570s");

    return 0;
}
