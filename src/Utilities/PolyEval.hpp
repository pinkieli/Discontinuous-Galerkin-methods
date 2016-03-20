#include <lapacke.h>

#ifndef PolyEval_HPP
#define PolyEval_HPP

double polyEval(double *GivenPoly, unsigned deg, double x)
{
	double Result = GivenPoly[deg];
	for(int i = deg-1;i>=0;i--)
		Result = Result*x+GivenPoly[i];
	return Result ;
}

#endif
