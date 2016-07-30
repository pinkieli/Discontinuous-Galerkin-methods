#include <lapacke.h>


#ifndef Ones_HPP
#define Ones_HPP

void ones(double *x, unsigned n)
{
    unsigned i;
    for(i=0; i<n; i++)
        x[i]    =   1.0;
}

void ones(double *A, unsigned m, unsigned n)
{
    unsigned i;
    for(i=0;i<m*n;i++)
            A[i] = 1.0;
}

#endif
