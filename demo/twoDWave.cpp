#include "src/includes.hpp"
#include <cmath>

#define U0  1.0

double U(double x, double y)
{
    return U0;
}

double V(double x, double y)
{
    return 0.0;
}

double initialConditions(double x, double y)
{
    return (exp(-8*(x*x+y*y)));
}

int main()
{
    unsigned    Nex =   10;
    unsigned    Ney =   10;
    unsigned    N   =   4;
    double L_start  =   -1;
    double L_end    =   1;
    double H_start  =   -1;
    double H_end    =   1;
    Field q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setVelocity(U,V);
    q.setInitialConditions(initialConditions);
    q.setSolver(5e-3,400);
    q.solve();
    q.plotSolution("time=400");

    return 0;
}
