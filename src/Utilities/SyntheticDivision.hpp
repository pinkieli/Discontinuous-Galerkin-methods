#include <lapacke.h>
#include <cstring>

#ifndef SyntheticDivision_HPP
#define SyntheticDivision_HPP

void syntheticDivision(double *GivenPoly, unsigned deg, double root)
{
    int i;
    double *Result;
    Result = new double [deg];


	Result[deg-1] = GivenPoly[deg];
	for(i = deg-2;i>=0;i--)
		Result[i] = Result[i+1]*root + GivenPoly[i+1];

    memcpy(GivenPoly,Result,deg*(sizeof(double)));
    GivenPoly[deg] = 0;

    delete [] Result;
    return ;
}


#endif
