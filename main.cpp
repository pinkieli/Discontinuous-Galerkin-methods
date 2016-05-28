#include "src/includes.hpp"
#include <cmath>
#include <string>
using namespace std;


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
    return (exp(-16*(x*x+y*y)));
}

double Depth(double x, double y)
{
    return 2.0;
}
int main()
{
    string name;
    unsigned    Nex =   20;
    unsigned    Ney =   20;
    unsigned    N   =   4 ;
    double L_start  =   -2;
    double L_end    =   2;
    double H_start  =   -2;
    double H_end    =   2;
    unsigned NTimeSteps =200;
    double dt  =   1e-5;
    unsigned i=0;


    ShallowWater q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setDepth(Depth);
    q.setInitialConditions(eta,U,V);
    q.setSolver(dt,NTimeSteps);
    for(i=0;i<50;i++)
    {
        name =  "t="+to_string(i*NTimeSteps*dt)+"s";
        q.plotSolution(0.0,3.0,name);
        q.solve();
    }
    name =  "t="+to_string(i*NTimeSteps*dt)+"s";
    q.plotSolution(0.0,3.0,name);

    return 0;
}
