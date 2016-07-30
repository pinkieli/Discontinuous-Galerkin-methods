#include "src/includes.hpp"


int main()
{
    unsigned N=2;
    double derivative[N+1][N+1];
    derivativeMatrix(*derivative,N);
    display(*derivative,N+1,N+1);

    return 0;
}
