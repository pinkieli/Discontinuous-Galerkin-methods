#include "../src/includes.hpp"
#include <cmath>
#include <string>
using namespace std;

double f(double x) {
    return cos(x);
}

int main()  {
    unsigned N = 1000;
    double L_start  =   -3.14159;
    double L_end    =    3.14159;
    double h        =   (L_end-L_start)/N;
    unsigned i;
    double X[N],Y[N];
    X[0]    =   L_start;
    Y[0]    =   f(X[0]);
    for(i=1;i<N;i++) {
        X[i]    =   X[i-1]  +   h;
        Y[i]    =   f(X[i]);
    }

    plot(X,Y,N,"Try","Density","CoSine");

    return 0;
}
