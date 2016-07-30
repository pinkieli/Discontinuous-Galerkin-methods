#include <lapacke.h>
#include <functional>
#define epsi 0.001
using namespace std;

#ifndef NewtonRaphson_H
#define NewtonRaphson_H

double newtonRaphson(function<double(double)> func ,double x)
{
	double deriv;
	double x_curr= x;
	double x_pre = x;
	do
	{
		x_pre = x_curr;
		deriv = (func(x_pre+epsi)-func(x_pre-epsi))/(2*epsi);
		x_curr = x_pre - (func(x_pre))/(deriv);
	}while(abs(x_curr - x_pre) >= 1e-6);
	return x_curr ;
}

#endif
