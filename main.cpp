#include "src/includes.hpp"
#include <cmath>

#define NEL 28

double U(double x, double y)
{
    return y;
}

double V(double x, double y)
{
    return -x;
}

double initialConditions(double x, double y)
{
    return (exp(-32*((x+0.5)*(x+0.5)+y*y)));;
}

double exactSolution(double x, double y)
{
    double x0   =   -0.5*cos(6.250);
    double y0   =   0.5*sin(6.250);;
    return (exp(-32*((x-x0)*(x-x0)+(y-y0)*(y-y0))));;
}

int main()
{
    unsigned    Nex =   NEL;
    unsigned    Ney =   NEL;
    unsigned    N   =   2;
    double L_start  =   -1;
    double L_end    =   1;
    double H_start  =   -1;
    double H_end    =   1;
    Field q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setVelocity(U,V);
    q.setInitialConditions(initialConditions);
    q.setSolver(0.0625*0.5*(1.0/NEL),100*2*(NEL));
    q.solve();
    q.plotSolution(0,1,"t=2sec");
    printf("%6.6f\tfor N = %d N_p = %d\n",q.L2Error(exactSolution),N,(Nex*Ney*(N+1)*(N+1)));

    return 0;
}
