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
}

double Depth(double x, double y)
{
    return 2.0;
}
int main()
{
    string name;
    unsigned    Nex =   20;
    unsigned    Ney =   5;
    unsigned    N   =   16 ;
    double L_start  =   -2;
    double L_end    =   2;
    double H_start  =   0;
    double H_end    =   1;
    unsigned NTimeSteps =50;
    double dt  =  1e-4;
    unsigned i;


    ShallowWater q(Nex,Ney,N);
    q.setDomain(L_start,L_end,H_start,H_end);
    q.setDepth(Depth);
    q.setInitialConditions(eta,U,V);
    q.setSolver(dt,NTimeSteps);
    for(i=0;i<2000;i++)
    {
        name =  "DamBreakN=16/DamBreak"+to_string(i)+".vtk";
        q.writeVTK(name);
        q.solve();
    }
    name =  "DamBreakN=16/DamBreak"+to_string(i)+".vtk";
    q.writeVTK(name);

    return 0;
}
