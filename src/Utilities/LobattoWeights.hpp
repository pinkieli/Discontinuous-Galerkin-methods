#include <lapacke.h>
#include <functional>

#include "LobattoNodes.hpp"
#include "PolyEval.hpp"
#include "PolyDeriv.hpp"
#include "LegendrePolynomial.hpp"

#ifndef LobattoWeights_HPP
#define LobattoWeights_HPP

void lobattoWeights(double *Weights, unsigned N)
{
    double *Poly, *Nodes;
    Poly    =   new double[N];
    Nodes   =   new double[N];
	legendrePolynomial(Poly,N-1);
    lobattoNodes(Nodes,N);

    function<double(double)> Eval;
    Weights[0]  = 2.0/((N)*(N-1));

    Eval = [&Poly,&N](double x){ return(polyEval(Poly,N-1,x));};

	for(int i=1;i<N-1;i++)
		Weights[i]    =   2/((N*(N-1))*(Eval(Nodes[i]))*(Eval(Nodes[i])));
	Weights[N-1] = (2.0/((N)*(N-1)));
}

#endif
