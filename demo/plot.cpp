#include "src/includes.hpp"
#include <cmath>

int main()
{
    int N =100;
    double X[N][N],Y[N][N],Z[N][N],x,y;
    double h = 2.0/N;
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            x   =   X[i][j]=-1  +   h*j;
            y   =   Y[i][j]=-1  +   h*i;
            Z[i][j] =       exp(-8*(x*x + y*y));
        }
    }
    plot(*X,*Y,*Z,N,N,"danda");

    return 0;
}
