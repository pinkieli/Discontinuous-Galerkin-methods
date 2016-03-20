#include "src/includes.hpp"

int main()
{
    double poly[4];
    poly[0]=-9;
    poly[1]=13.5;
    poly[2]=-6.5;
    poly[3]=1;
    syntheticDivision(poly,4,1.5); 
    display(poly,3);




    return 0;
}
