#include "src/includes.hpp"


int main()
{
    unsigned N = 1 ;
    double DerivativeMatrixX[(N+1)*(N+1)][(N+1)*(N+1)];
    twoDDerivativeMatrixX(*DerivativeMatrixX,N); 
    display(*DerivativeMatrixX,(N+1)*(N+1),(N+1)*(N+1));

    printf("\n================================\n\n");

    double DerivativeMatrixY[(N+1)*(N+1)][(N+1)*(N+1)];
    twoDDerivativeMatrixY(*DerivativeMatrixY,N); 
    display(*DerivativeMatrixY,(N+1)*(N+1),(N+1)*(N+1));
    return 0;
}
