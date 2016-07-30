#include "src/includes.hpp"


double check(double x)
{
    return (cos(x));
}


int main()
{
    unsigned N = 4;
    double integral;
    integral    =   lobattoIntegration(0,0.5*3.14159,N,check);
    printf("%6.2f\n",integral);
    return 0;
}
