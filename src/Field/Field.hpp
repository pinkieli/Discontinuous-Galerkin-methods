#include <lapacke.h>
#include <functional>
#include <cstring>
#include <cblas.h>
#include "../Utilities/LobattoNodes.hpp"
#include "../Utilities/Inverse.hpp"
#include "../Elemental/MassMatrix.hpp"
#include "../Elemental/DerivativeMatrix.hpp"
using namespace std;

#ifndef Field_HPP
#define Field_HPP

class Field
{
private:
    unsigned Nex,Ney,N;
    double L_start,L_end,H_start,H_end;

    double ***ConsVariable,***Rate;
    double ***X, ***Y;
    double ***U, ***V;

    double dt;
    unsigned NTimeSteps;

public:

    Field(unsigned , unsigned, unsigned );

    void setDomain(double , double , double , double );
    void setVelocity(function<double(double,double)> , function<double(double,double)> );
    void setInitialConditions(function<double(double,double)> );
    void setSolver(double , unsigned );

    void solve();

};

Field::Field(unsigned Nx, unsigned Ny, unsigned n)
{
    unsigned i,j;
    Nex =   Nx;
    Ney =   Ny;
    N   =   n;
    ConsVariable    =   new double** [Ney];
    for( i = 0; i< Ney; i++)
    {
        ConsVariable[i]    =   new double* [Nex];
        for(j=0;j<Nex;j++)
            ConsVariable[i][j]  =   new double[(N+1)*(N+1)];
    }

    Rate            =   new double** [Ney];
    for( i = 0; i< Ney; i++)
    {
        Rate[i]     =   new double* [Nex];
        for(j=0;j<Nex;j++)
            Rate[i][j]      =   new double[(N+1)*(N+1)];
    }

    X    =   new double** [Ney];
    for( i = 0; i< Ney; i++)
    {
        X[i]    =   new double* [Nex];
        for(j=0;j<Nex;j++)
            X[i][j]  =   new double[(N+1)*(N+1)];
    }

    Y    =   new double** [Ney];
    for( i = 0; i< Ney; i++)
    {
        Y[i]    =   new double* [Nex];
        for(j=0;j<Nex;j++)
            Y[i][j]  =   new double[(N+1)*(N+1)];
    }

}

void Field::setDomain(double L1, double L2, double L3, double L4)
{
    unsigned i,j,k1,k2;
    L_start =   L1;
    L_end   =   L2;
    H_start =   L3;
    H_end   =   L4;

    double Xstart,Ystart;
    double Xend,Yend;
    double dx,dy;
    double Nodes[N+1];
    lobattoNodes(Nodes,N+1);


    dx  =   (L_end-L_start)/Nex;
    dy  =   (H_end-H_start)/Ney;

    Ystart  =   H_start;
    Yend    =   Ystart+dy;

    for( i = 0;i<Ney; i++)
    {
        Xstart  =   L_start;
        Xend    =   L_end;
        for( j=0 ; j<Nex; j++)
        {
            for(k1=0;k1<=N;k1++)
            {
                for(k2=0;k2<=N;k2++)
                {
                    X[i][j][k1*(N+1)+k2]=Xstart + 0.5*dx*Nodes[k2];
                    Y[i][j][k1*(N+1)+k2]=Ystart + 0.5*dy*Nodes[k1];
                }
            }
            Xstart+=dx;
            Xend+=dx;
        }
        Ystart+=dy;
        Yend+=dy;
    }
    return ;
}

void Field::setVelocity(function<double(double,double)> A, function<double(double,double)> B)
{
    unsigned i,j,k;

    for ( i=0;i<Ney;i++)
    {
        for ( j=0; j<Nex; j++)
        {
            for (k=0;k<((N+1)*(N+1));k++)
            {
                U[i][j][k]  =   A(X[i][j][k],Y[i][j][k]);
                V[i][j][k]  =   B(X[i][j][k],Y[i][j][k]);
            }
        }
    }

    return ;
}

void Field::setInitialConditions(function<double(double,double)> A)
{
    unsigned i,j,k;

    for ( i=0;i<Ney;i++)
        for ( j=0; j<Nex; j++)
            for (k=0;k<((N+1)*(N+1));k++)
                ConsVariable[i][j][k]  =   A(X[i][j][k],Y[i][j][k]);
    return ;
}

void Field::setSolver(double a, unsigned b)
{
    dt          =   a;
    NTimeSteps  =   b;
}




void Field::solve()
{
    unsigned i,j,k;
    double MassMatrix[(N+1)*(N+1)][(N+1)*(N+1)];
    double MassInverse[(N+1)*(N+1)][(N+1)*(N+1)];
    double DerivativeMatrixX[(N+1)*(N+1)][(N+1)*(N+1)], DerivativeMatrixY[(N+1)*(N+1)][(N+1)*(N+1)];

    twoDMassMatrix(*MassMatrix,N);
    inverse(*MassMatrix,*MassInverse,(N+1)*(N+1));

    twoDDerivativeMatrixX(*DerivativeMatrixX,N);
    twoDDerivativeMatrixY(*DerivativeMatrixY,N);

    double ***RHS, ***PreviousValue, ***f_preX,***f_preY;
    RHS             =   new double** [Ney];
    PreviousValue   =   new double** [Ney];
    f_preX          =   new double** [Ney];
    f_preY          =   new double** [Ney];
    for( i = 0; i< Ney; i++)
    {
        f_preX[i]           =   new double* [Nex];
        f_preY[i]           =   new double* [Nex];
        RHS[i]              =   new double* [Nex];
        PreviousValue[i]    =   new double* [Nex];
        for(j=0;j<Nex;j++)
        {
            f_preX[i][j]        =   new double[(N+1)*(N+1)];
            f_preY[i][j]        =   new double[(N+1)*(N+1)];
            RHS[i][j]           =   new double[(N+1)*(N+1)];
            PreviousValue[i][j] =   new double[(N+1)*(N+1)];
        }
    }

    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            memcpy(RHS[i][j],ConsVariable[i][j],(N+1)*(N+1)*sizeof(double));
            memcpy(PreviousValue[i][j],ConsVariable[i][j],(N+1)*(N+1)*sizeof(double));
            for(k=0;k<(N+1)*(N+1);k++)
            {
                f_preX[i][j][k]  = PreviousValue[i][j][k]*U[i][j][k];
                f_preY[i][j][k]  = PreviousValue[i][j][k]*V[i][j][k];
            }
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),1,*DerivativeMatrixX,(N+1)*(N+1),f_preX[i][j],1,0,RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),1,*DerivativeMatrixY,(N+1)*(N+1),f_preY[i][j],1,1,RHS[i][j],1);
        }
    }

    /*HERE FLUX WILL DO THE CHANGES TO THE RHS USING THE PREVIOUS VALUE AND RHS*/
    /*What I am planning right now is continue to build up the flux matrices as in the previous case and then figure out a way to deal with the non - contiguity.*/
    /*OR there could be a way of iterating through the sides and getting the result like done by Giraldo.*/

    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),1,*MassInverse,(N+1)*(N+1),RHS[i][j],1,0,Rate[i][j],1);

    return ;
}

#endif
