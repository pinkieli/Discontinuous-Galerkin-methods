#include "src/includes.hpp"

int main()
{
    unsigned N = 6;
    double Weights[N];
    lobattoWeights(Weights,N);
    display(Weights,N);
    return 0;
}
