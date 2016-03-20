#include "src/includes.hpp"

int main()
{
    double poly[4];
    poly[0]=-9;
    poly[1]=13.5;
    poly[2]=-6.5;
    poly[3]=1;
    double x =polyEval(poly,3,20.23); 
    printf("%6.2f\n",x);




    return 0;
}
