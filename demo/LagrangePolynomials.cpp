#include "src/includes.hpp"


int main()
{
    unsigned N=5;
    double Poly[N+1][N+1];
    lagrangePolynomials(*Poly,N);
    display(*Poly,N+1,N+1);
    return 0;
}
