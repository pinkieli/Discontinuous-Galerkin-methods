#include "src/includes.hpp"


int main()
{
    unsigned N = 2 ;
    double FluxMatrix[(N+1)][(N+1)];
    fluxMatrix(*FluxMatrix,N); 
    display(*FluxMatrix,(N+1),(N+1));

    return 0;
}
