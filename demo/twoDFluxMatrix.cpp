#include "src/includes.hpp"


int main()
{
    unsigned N = 1 ;
    double FluxMatrix[(N+1)*(N+1)][(N+1)*(N+1)];
    twoDFluxMatrix1(*FluxMatrix,N);
    display(*FluxMatrix,(N+1)*(N+1),(N+1)*(N+1));
    printf("\n=========================\n\n");

    twoDFluxMatrix2(*FluxMatrix,N);
    display(*FluxMatrix,(N+1)*(N+1),(N+1)*(N+1));
    printf("\n=========================\n\n");

    twoDFluxMatrix3(*FluxMatrix,N);
    display(*FluxMatrix,(N+1)*(N+1),(N+1)*(N+1));
    printf("\n=========================\n\n");

    twoDFluxMatrix4(*FluxMatrix,N);
    display(*FluxMatrix,(N+1)*(N+1),(N+1)*(N+1));

    return 0;
}
