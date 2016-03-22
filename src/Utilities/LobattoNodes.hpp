#include <lapacke.h>
#include "LegendrePolynomial.hpp"
#include "PolyEval.hpp"
#include "PolyDeriv.hpp"
#include "SyntheticDivision.hpp"
#include "NewtonRaphson.hpp"
#include "Display.hpp"
#include <algorithm>
#include <functional>
using namespace std;

#ifndef LobattoNodes_HPP
#define LobattoNodes_HPP

void lobattoNodes(double *Nodes, unsigned N)
{
    double *Poly,*DerivedPoly;
    Poly        =   new double[N];
    DerivedPoly =   new double[N-1];
    double root;
	double InitialGuess=-1.0;
	unsigned deg=N-2;///To store the degree at each state of the polynomial.
    legendrePolynomial(Poly,N-1);
    polyDeriv(Poly,DerivedPoly,N-1);

    function<double(double)> eval;
	Nodes[0]   =   -1;
	for(int i=1;i<=N-2;i++)
	{
		eval = [&DerivedPoly,&deg](double x){ return(polyEval(DerivedPoly,deg,x));};
		root = newtonRaphson(eval,InitialGuess);
		syntheticDivision(DerivedPoly,deg,root);
        --deg;
		Nodes[i]  =   root;
	}
	Nodes[N-1] =   1;
	sort(Nodes,Nodes+N);

    delete[] Poly;
    delete[] DerivedPoly;
}


#endif
