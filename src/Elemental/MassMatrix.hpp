#include <lapacke.h>
#include <functional>

using namespace std;

#ifndef MassMatrix_HPP
#define MassMatrix_HPP

void massMatrix(double *MassMatrix,unsigned N)
{
    vector< vector<double> > MassMatrix;
    vector< vector<double> > LagrangePolynomials = lagrangePolynomials(Points);
    unsigned i,j;///Counters for the loop.
    function<double(double)> eval;
    for(i=0;i<n;i++)
    {

        for(j=0;j<n;j++)
        {
            eval = [&LagrangePolynomials,&i,&j](double x){  return ((polyEval(LagrangePolynomials[i],x)*polyEval(LagrangePolynomials[j],x)));};
            MassMatrix[i][j] = lobattoIntegration(start,end,n,eval);
        }
    }
}

#endif
