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
    if(x<0.0)
        return 1.0;
    else
        return 0.0;
    //return (exp(-16*((x)*(x)+y*y)));
}

double Depth(double x, double y)
{
    return 2.0;
}
int main()
{
    string name;
    unsigned    Nex =   100;
    unsigned    Ney =   2;
    unsigned    N   =   8 ;
    double L_start  =   -2;
    double L_end    =   2;
    double H_start  =   0.0;
    double H_end    =   0.1;
    unsigned NTimeSteps =50;
    double dt  =  1e-4;
    unsigned i;


    ShallowWater q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setDepth(Depth);
    q.setInitialConditions(eta,U,V);
    q.setSolver(dt,NTimeSteps);
    for(i=0;i<1000;i++)
    {
        name =  "t="+to_string(i*NTimeSteps*dt)+"s";
        q.plotBoundary(1.5,3.5,name);
        q.solve();
    }
    name =  "t="+to_string(i*NTimeSteps*dt)+"s";
    q.plotBoundary(1.5,3.5,name);

    return 0;
}
