#include "src/includes.hpp"


int main()
{
    unsigned N = 1 ;
    double MassMatrix[(N+1)*(N+1)][(N+1)*(N+1)];
    twoDMassMatrix(*MassMatrix,N); 
    display(*MassMatrix,(N+1)*(N+1),(N+1)*(N+1));

    return 0;
}
