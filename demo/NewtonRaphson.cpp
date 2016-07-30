#include "src/includes.hpp"

int main()
{
    double poly[4];
    poly[0]=-9;
    poly[1]=13.5;
    poly[2]=-6.5;
    poly[3]=1;
    function<double(double)> eval;
    eval = [&poly](double x){return polyEval(poly,3,x);};
    double y =newtonRaphson(eval,0);    
    printf("%6.2f\n",y);





    return 0;
}
