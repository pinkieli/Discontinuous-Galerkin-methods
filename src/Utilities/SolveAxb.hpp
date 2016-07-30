#include <lapacke.h>
#include <cstring>
#include <cstdio>
#include "Display.hpp"

#ifndef SolveAxb_HPP
#define SolveAxb_HPP

void solveAxb(double *A, double *x, double *b, unsigned N)
{
    double *B    = new double[N*N];
    double *c    = new double[N];
    memcpy(B,A,N*N*sizeof(double));
    memcpy(c,b,N*sizeof(double));
    int ipiv[N];
    int info;

    info    =   LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,B,N,ipiv,c,1);
    if(info!=0)
        fprintf(stderr,"The Linear solve `Ax=b` was not succesful.\n");
   
    memcpy(x,c,N*sizeof(double));

    delete[] B;
    delete[] c;
    return ;
}



#endif
