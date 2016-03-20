#include "src/includes.hpp"

int main()
{
    double A[3][3];
    double Ainv[3][3];
    A[0][0]=1;
    A[0][1]=1;
    A[0][2]=7;

    A[1][0]=1;
    A[1][1]=2;
    A[1][2]=3;
    
    A[2][0]=2;
    A[2][1]=3;
    A[2][2]=5;

    inverse(*A,*Ainv,3);

        

    display(*Ainv,3,3);

    return 0;
}
