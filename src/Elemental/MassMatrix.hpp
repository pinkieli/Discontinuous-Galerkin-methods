#include <lapacke.h>
#include <functional>
#include <cstring>
#include "../Utilities/LagrangePolynomials.hpp"
#include "../Utilities/PolyEval.hpp"
#include "../Utilities/LobattoIntegration.hpp"

using namespace std;

#ifndef MassMatrix_HPP
#define MassMatrix_HPP

void massMatrix(double *MassMatrix,unsigned N)
{
    double Poly[N+1][N+1];
    double **poly;
    poly   =   new double*[N+1];
    lagrangePolynomials(*Poly,N);
    unsigned i,j;

    for(i=0;i<=N;i++)
    {
        poly[i] =   new double[N+1];
        memcpy(poly[i],Poly[i],(N+1)*sizeof(double));
    }

    function<double(double)> eval;
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++)
        {
            eval = [&poly,&i,&j,&N](double x){return (((polyEval(poly[i],N,x))*(polyEval(poly[j],N,x))));};
            MassMatrix[i*(N+1)+j] = lobattoIntegration(-1.0,1.0,N+1,eval);
        }
    }
    return ;
}
#endif
