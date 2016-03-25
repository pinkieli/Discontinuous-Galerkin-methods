#include <lapacke.h>
#include "../Utilities/LagrangePolynomials.hpp"

#ifndef FluxMatrix_HPP
#define FluxMatrix_HPP

void fluxMatrix(double *FluxMatrix, unsigned N)
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
            FluxMatrix[i*(N+1)+j] = eval(1)-eval(-1);
        }
    }

    for(i=0;i<=N;i++)
    {
        delete[] poly[i];
    }

    delete[] poly;

    return ;
}


#endif
