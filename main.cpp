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
    if((x<0.2)&&(x>-0.2)&&(y<0.2)&&(y>-0.2))
        return 1;
    else
        return 0;
}

int main()
{
    unsigned    Nex =   20;
    unsigned    Ney =   20;
    unsigned    N   =   8;
    double L_start  =   -1;
    double L_end    =   1;
    double H_start  =   -1;
    double H_end    =   1;
    Field q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setVelocity(U,V);
    q.setInitialConditions(initialConditions);
    q.setSolver(5e-4,0);
    q.solve();
    q.plotSolution(0,1,"t=0sec");

    return 0;
}
