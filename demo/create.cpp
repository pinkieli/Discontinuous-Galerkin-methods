#include "src/includes.hpp"


int main()
{
    double ***A;
    A   =   create3D(2,2,2);
    unsigned i,j,k;
    for(i=0;i<2;i++)
        for(j=0;j<2;j++)
            for(k=0;k<2;k++)
                A[i][j][k]  =   i+j+k;

    i=0;
    for(j=0;j<2;j++)
    {
        for(k=0;k<2;k++)
            printf("%6.2f\t",A[i][j][k] );
        printf("\n");
    }
    printf("\n\n\n");
    i=1;
    for(j=0;j<2;j++)
    {
        for(k=0;k<2;k++)
            printf("%6.2f\t",A[i][j][k] );
        printf("\n");
    }
    return 0;
}
