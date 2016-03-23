#include "src/includes.hpp"


int main()
{
    unsigned N=2;
    double Mass[N+1][N+1];
    massMatrix(*Mass,N);
    display(*Mass,N+1,N+1);

    return 0;
}
