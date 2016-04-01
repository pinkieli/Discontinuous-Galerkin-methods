#ifndef Create_HPP
#define Create_HPP

double*** create3D(unsigned X,unsigned  Y,unsigned  Z)
{
    double ***A;
    unsigned i,j;
    A    =   new double** [X];
    for( i = 0; i< X; i++)
    {
        A[i]    =   new double* [Y];
        for(j=0;j<Y;j++)
            A[i][j]  =   new double[Z];
    }
    return A;
}

#endif
