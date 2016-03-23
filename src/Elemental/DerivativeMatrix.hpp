#include <lapacke.h>
#include <functional>
#include <cstring>
#include "../Utilities/LagrangePolynomials.hpp"
#include "../Utilities/PolyEval.hpp"
#include "../Utilities/PolyDeriv.hpp"
#include "../Utilities/LobattoIntegration.hpp"

using namespace std;

#ifndef DerivativeMatrix_HPP
#define DerivativeMatrix_HPP

void derivativeMatrix(double *DerivativeMatrix,unsigned N)
{
    double Poly[N+1][N+1];
    double **poly, **deriv;
    poly   =   new double*[N+1];
    deriv  =   new double*[N+1];
    lagrangePolynomials(*Poly,N);
    unsigned i,j;

    for(i=0;i<=N;i++)
    {
        poly[i] =   new double[N+1];
        deriv[i]=   new double[N];
        memcpy(poly[i],Poly[i],(N+1)*sizeof(double));
        polyDeriv(poly[i],deriv[i],N+1);
    }

    function<double(double)> eval;
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++)
        {
            eval = [&poly,&deriv,&i,&j,&N](double x){return (((polyEval(poly[i],N,x))*(polyEval(deriv[j],N-1,x))));};
            DerivativeMatrix[i*(N+1)+j] = lobattoIntegration(-1.0,1.0,N+1,eval);
        }
    }

    for(i=0;i<=N;i++)
    {
        delete[] poly[i];
        delete[] deriv[i];
    }


    delete[] poly;
    delete[] deriv;
    return ;
}
#endif
