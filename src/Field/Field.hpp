#include <lapacke.h>
#include <functional>
#include <cstring>
#include <cblas.h>
#include <cmath>
#include <string>
#include "../Utilities/LobattoNodes.hpp"
#include "../Utilities/Inverse.hpp"
#include "../Utilities/Zeros.hpp"
#include "../Utilities/Create.hpp"
#include "../Utilities/Plot.hpp"
#include "../Elemental/MassMatrix.hpp"
#include "../Elemental/DerivativeMatrix.hpp"
#include "../Elemental/FluxMatrix.hpp"

#define MAX(a, b)(a>b?a:b)
#define ABS(a)(a>0?a:(-a))

using namespace std;

#ifndef Field_HPP
#define Field_HPP



class Field
{
private:
    /**
    * Nex   : This denotes the number of elements in the X- direction.
    * Ney   : This denotes the number of elements in the Y- direction.
    * N     : This variable denotes the order of the polynomial.
    * Note  : This code considers same order of polynomail approximation in both X- and Y- directions.
    */
    unsigned Nex,Ney,N;

    /**
    * The variables listed below are for defining the domain of the problem.
    * L_start: This denotes the X- co-ordinate of the starting point..
    * L_end:   This denotes the X-coordinate of the ending point
    * H_start: This denotes the Y- co-ordinate of the starting point..
    * H_end:   This denotes the Y-coordinate of the ending point
    */

    double L_start,L_end,H_start,H_end;
    /**
    * double ***ConsVariable: This variable is a dynamically allocated 3-D array for the Conservation variable on which the Hyperbolic equation is acting.
    */

    double ***ConsVariable,***Rate,***RHS,***XFlux,***YFlux;
    double ***X, ***Y;
    double ***U, ***V;

    double dt;
    unsigned NTimeSteps;

    double dx,dy;

    double Lambda;

public:

    Field(unsigned , unsigned, unsigned );

    void setDomain(double , double , double , double );
    void setVelocity(function<double(double,double)> , function<double(double,double)> );
    void setInitialConditions(function<double(double,double)> );
    void setSolver(double , unsigned );
    void computeLambda();
    void computeFlux(double*** , double*** );
    void operateDerivative(double ***, double ***, double *, double * );
    void computeNumericalFlux(double ***, double ***);
    void operateFlux(double *,double *,double *,double *);
    void operateInvereseMass(double *);
    void copyField(double ***);
    void plotSolution(string );
    void solve();
};

Field::Field(unsigned Nx, unsigned Ny, unsigned n)
{
    Nex =   Nx;
    Ney =   Ny;
    N   =   n;
    ConsVariable    =   create3D(Ney,Nex,(N+1)*(N+1));
    Rate            =   create3D(Ney,Nex,(N+1)*(N+1));
    RHS             =   create3D(Ney,Nex,(N+1)*(N+1));
    XFlux           =   create3D(Ney,Nex,(N+1)*(N+1));
    YFlux           =   create3D(Ney,Nex,(N+1)*(N+1));
    X               =   create3D(Ney,Nex,(N+1)*(N+1));
    Y               =   create3D(Ney,Nex,(N+1)*(N+1));
    U               =   create3D(Ney,Nex,(N+1)*(N+1));
    V               =   create3D(Ney,Nex,(N+1)*(N+1));
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
    double Nodes[N+1];
    lobattoNodes(Nodes,N+1);

    dx  =   (L_end-L_start)/Nex;
    dy  =   (H_end-H_start)/Ney;

    Ystart  =   H_start;
    Yend    =   Ystart+dy;

    for( i = 0;i<Ney; i++)
    {
        Xstart  =   L_start;
        Xend    =   Xstart + dx;
        for( j=0 ; j<Nex; j++)
        {
            for(k1=0;k1<=N;k1++)
            {
                for(k2=0;k2<=N;k2++)
                {
                    X[i][j][k1*(N+1)+k2]=0.5*(Xstart+Xend) + 0.5*dx*Nodes[k2];
                    Y[i][j][k1*(N+1)+k2]=0.5*(Ystart+Yend) + 0.5*dy*Nodes[k1];
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

void Field::computeLambda()
{

    Lambda = 0;
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                Lambda  =   MAX(Lambda,ABS(V[i][j][k]));
                Lambda  =   MAX(Lambda,ABS(U[i][j][k*(N+1)]));
            }
        }
    }

    j   =   Nex-1;

    for(i=0;i<Ney;i++)
        for(k=0;k<=N;k++)
            Lambda  =   MAX(Lambda,ABS(U[i][j][k*(N+1)+N]));

    i   =   Ney-1;

    for(j=0;j<Nex;j++)
        for(k=0;k<=N;k++)
            Lambda  =   MAX(Lambda,ABS(V[i][j][k + N*(N+1) ]));



    return ;
}

void Field::computeFlux(double ***f_preX, double ***f_preY )
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<((N+1)*(N+1));k++)
            {
                f_preX[i][j][k]  = ConsVariable[i][j][k]*U[i][j][k];
                f_preY[i][j][k]  = ConsVariable[i][j][k]*V[i][j][k];
            }
        }
    }
    return ;
}

void Field::operateDerivative(double ***f_preX, double ***f_preY , double *DerivativeMatrixX, double *DerivativeMatrixY )
{
    unsigned i,j;

    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),f_preX[i][j],1,0,RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),f_preY[i][j],1,1,RHS[i][j],1);
        }
    }

    return ;
}

void Field::computeNumericalFlux(double ***f_preX, double ***f_preY )
{
    unsigned i,j,k;

    for(i=0;i<Ney;i++)
    {
        for(j=1;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                XFlux[i][j][k*(N+1)]    =   0.5*(f_preX[i][j][k*(N+1)] + f_preX[i][j-1][k*(N+1) + N] - Lambda*(ConsVariable[i][j][k*(N+1)]-ConsVariable[i][j-1][k*(N+1)+N]) );
                XFlux[i][j-1][k*(N+1)+N]=   0.5*(f_preX[i][j][k*(N+1)] + f_preX[i][j-1][k*(N+1) + N] - Lambda*(ConsVariable[i][j][k*(N+1)]-ConsVariable[i][j-1][k*(N+1)+N]) );
            }
        }
    }

    j=0;

    for(i=0;i<Ney;i++)
    {
        for(k=0;k<=N;k++)
        {
            XFlux[i][j][k*(N+1)]        =   0.5*(f_preX[i][j][k*(N+1)] + f_preX[i][Nex-1][k*(N+1) + N] - Lambda*(ConsVariable[i][j][k*(N+1)]-ConsVariable[i][Nex-1][k*(N+1)+N]) );
            XFlux[i][Nex-1][k*(N+1)+N]  =   0.5*(f_preX[i][j][k*(N+1)] + f_preX[i][Nex-1][k*(N+1) + N] - Lambda*(ConsVariable[i][j][k*(N+1)]-ConsVariable[i][Nex-1][k*(N+1)+N]) );
        }
    }

    for(i=1;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                YFlux[i][j][k]              =   0.5*(f_preY[i][j][k] + f_preY[i-1][j][k + N*(N+1)] - Lambda*(ConsVariable[i][j][k] - ConsVariable[i-1][j][k + N*(N+1)]) );
                YFlux[i-1][j][k + N*(N+1)]  =   0.5*(f_preY[i][j][k] + f_preY[i-1][j][k + N*(N+1)] - Lambda*(ConsVariable[i][j][k] - ConsVariable[i-1][j][k + N*(N+1)]) );
            }
        }
    }

    i=0;

    for(j=0;j<Nex;j++)
    {
        for(k=0;k<=N;k++)
        {
            YFlux[i][j][k]                  =   0.5*(f_preY[i][j][k] + f_preY[Ney-1][j][k + N*(N+1)] - Lambda*(ConsVariable[i][j][k] - ConsVariable[Ney-1][j][k + N*(N+1)]) );
            YFlux[Ney-1][j][k + N*(N+1)]    =   0.5*(f_preY[i][j][k] + f_preY[Ney-1][j][k + N*(N+1)] - Lambda*(ConsVariable[i][j][k] - ConsVariable[Ney-1][j][k + N*(N+1)]) );
        }
    }

    return ;
}

void Field::operateFlux(double *Flux1,double *Flux2,double *Flux3,double *Flux4)
{
    unsigned i,j;

    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dy,Flux2,(N+1)*(N+1),XFlux[i][j],1,1,RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dy,Flux4,(N+1)*(N+1),XFlux[i][j],1,1,RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dx,Flux3,(N+1)*(N+1),YFlux[i][j],1,1,RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dx,Flux1,(N+1)*(N+1),YFlux[i][j],1,1,RHS[i][j],1);
        }
    }

    return ;
}

void Field::operateInvereseMass(double *MassInverse)
{
    unsigned i,j;
    double alpha    =   4/(dx*dy);
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),alpha,MassInverse,(N+1)*(N+1),RHS[i][j],1,0,Rate[i][j],1);

    return ;
}

void Field::copyField(double ***q)
{
    unsigned i,j;
    for( i = 0; i< Ney; i++)
        for(j=0;j<Nex;j++)
            memcpy(q[i][j],ConsVariable[i][j],(N+1)*(N+1)*sizeof(double));
    return ;
}

void Field::plotSolution(string s)
{
    double CG[Ney*N+1][Nex*N+1],CGX[Ney*N+1][Nex*N+1],CGY[Ney*N+1][Nex*N+1];
    zeros(*CG,Ney*N+1,Nex*N+1);
    zeros(*CGX,Ney*N+1,Nex*N+1);
    zeros(*CGY,Ney*N+1,Nex*N+1);
    unsigned i,j,k1,k2;
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k1=1;k1<N;k1++)
            {
                k2=0;
                CG[i*N+k1][j*N+k2]+=0.5*ConsVariable[i][j][k1*(N+1)+k2];

                for(k2=1;k2<N;k2++)
                    CG[i*N+k1][j*N+k2]=ConsVariable[i][j][k1*(N+1)+k2];

                k2=N;
                CG[i*N+k1][j*N+k2]+=0.5*ConsVariable[i][j][k1*(N+1)+k2];
            }

            k1=0;
            k2=0;
            CG[i*N+k1][j*N+k2]+=0.25*ConsVariable[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CG[i*N+k1][j*N+k2]+=0.5*ConsVariable[i][j][k1*(N+1)+k2];
            k2=N;
            CG[i*N+k1][j*N+k2]+=0.25*ConsVariable[i][j][k1*(N+1)+k2];

            k1=N;
            k2=0;
            CG[i*N+k1][j*N+k2]+=0.25*ConsVariable[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CG[i*N+k1][j*N+k2]+=0.5*ConsVariable[i][j][k1*(N+1)+k2];
            k2=N;
            CG[i*N+k1][j*N+k2]+=0.25*ConsVariable[i][j][k1*(N+1)+k2];
        }
    }

    for(i=0;i<(Ney*N+1);i++)
    {
        CG[i][0]    =2*CG[i][0];
        CG[i][Nex*N]=2*CG[i][Nex*N];
    }

    for(j=0;j<(Nex*N+1);j++)
    {
        CG[0][j]    =2*CG[0][j];
        CG[Ney*N][j]=2*CG[Ney*N][j];
    }

    /*Starting the CG form for X*/
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k1=1;k1<N;k1++)
            {
                k2=0;
                CGX[i*N+k1][j*N+k2]+=0.5*X[i][j][k1*(N+1)+k2];

                for(k2=1;k2<N;k2++)
                    CGX[i*N+k1][j*N+k2]=X[i][j][k1*(N+1)+k2];

                k2=N;
                CGX[i*N+k1][j*N+k2]+=0.5*X[i][j][k1*(N+1)+k2];
            }

            k1=0;
            k2=0;
            CGX[i*N+k1][j*N+k2]+=0.25*X[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CGX[i*N+k1][j*N+k2]+=0.5*X[i][j][k1*(N+1)+k2];
            k2=N;
            CGX[i*N+k1][j*N+k2]+=0.25*X[i][j][k1*(N+1)+k2];

            k1=N;
            k2=0;
            CGX[i*N+k1][j*N+k2]+=0.25*X[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CGX[i*N+k1][j*N+k2]+=0.5*X[i][j][k1*(N+1)+k2];
            k2=N;
            CGX[i*N+k1][j*N+k2]+=0.25*X[i][j][k1*(N+1)+k2];
        }
    }

    for(i=0;i<(Ney*N+1);i++)
    {
        CGX[i][0]    =2*CGX[i][0];
        CGX[i][Nex*N]=2*CGX[i][Nex*N];
    }

    for(j=0;j<(Nex*N+1);j++)
    {
        CGX[0][j]    =2*CGX[0][j];
        CGX[Ney*N][j]=2*CGX[Ney*N][j];
    }
    /*Starting the CG form for Y*/
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k1=1;k1<N;k1++)
            {
                k2=0;
                CGY[i*N+k1][j*N+k2]+=0.5*Y[i][j][k1*(N+1)+k2];

                for(k2=1;k2<N;k2++)
                    CGY[i*N+k1][j*N+k2]=Y[i][j][k1*(N+1)+k2];

                k2=N;
                CGY[i*N+k1][j*N+k2]+=0.5*Y[i][j][k1*(N+1)+k2];
            }

            k1=0;
            k2=0;
            CGY[i*N+k1][j*N+k2]+=0.25*Y[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CGY[i*N+k1][j*N+k2]+=0.5*Y[i][j][k1*(N+1)+k2];
            k2=N;
            CGY[i*N+k1][j*N+k2]+=0.25*Y[i][j][k1*(N+1)+k2];

            k1=N;
            k2=0;
            CGY[i*N+k1][j*N+k2]+=0.25*Y[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CGY[i*N+k1][j*N+k2]+=0.5*Y[i][j][k1*(N+1)+k2];
            k2=N;
            CGY[i*N+k1][j*N+k2]+=0.25*Y[i][j][k1*(N+1)+k2];
        }
    }

    for(i=0;i<(Ney*N+1);i++)
    {
        CGY[i][0]    =2*CGY[i][0];
        CGY[i][Nex*N]=2*CGY[i][Nex*N];
    }

    for(j=0;j<(Nex*N+1);j++)
    {
        CGY[0][j]    =2*CGY[0][j];
        CGY[Ney*N][j]=2*CGY[Ney*N][j];
    }

    plot(*CGX,*CGY,*CG,Ney*N+1,Nex*N+1,s);

    return ;
}

void Field::solve()
{
    unsigned i,j,t;
    double MassMatrix[(N+1)*(N+1)][(N+1)*(N+1)];
    double MassInverse[(N+1)*(N+1)][(N+1)*(N+1)];
    double Flux1[(N+1)*(N+1)][(N+1)*(N+1)],Flux2[(N+1)*(N+1)][(N+1)*(N+1)],Flux3[(N+1)*(N+1)][(N+1)*(N+1)],Flux4[(N+1)*(N+1)][(N+1)*(N+1)];
    double DerivativeMatrixX[(N+1)*(N+1)][(N+1)*(N+1)], DerivativeMatrixY[(N+1)*(N+1)][(N+1)*(N+1)];

    twoDMassMatrix(*MassMatrix,N);
    inverse(*MassMatrix,*MassInverse,(N+1)*(N+1));

    twoDDerivativeMatrixX(*DerivativeMatrixX,N);
    twoDDerivativeMatrixY(*DerivativeMatrixY,N);

    twoDFluxMatrix1(*Flux1,N);
    twoDFluxMatrix2(*Flux2,N);
    twoDFluxMatrix3(*Flux3,N);
    twoDFluxMatrix4(*Flux4,N);

    double ***f_preX,***f_preY,***q0,***q1,***q2;
    f_preX  =   create3D(Ney,Nex,((N+1)*(N+1)));
    f_preY  =   create3D(Ney,Nex,((N+1)*(N+1)));
    q0      =   create3D(Ney,Nex,((N+1)*(N+1)));
    q1      =   create3D(Ney,Nex,((N+1)*(N+1)));
    q2      =   create3D(Ney,Nex,((N+1)*(N+1)));

    computeLambda();

    for(t=0;t<NTimeSteps;t++)
    {
        printf("Time t=%6.3f\tCourant Number=%6.2f\n",(t+1.0)*dt,Lambda*dt*(1.0/dy+1.0/dx) );
        copyField(q0);
        computeFlux(f_preX,f_preY);
        computeNumericalFlux(f_preX, f_preY);
        operateDerivative(f_preX,f_preY,*DerivativeMatrixX,*DerivativeMatrixY);
        operateFlux(*Flux1,*Flux2,*Flux3,*Flux4);
        operateInvereseMass(*MassInverse);
        //q (1) = q (0) + ∆tR(q (0));
        for(i=0;i<Ney;i++)
            for(j=0;j<Nex;j++)
                cblas_daxpy((N+1)*(N+1),dt,Rate[i][j],1,ConsVariable[i][j],1);


        copyField(q1);
        computeFlux(f_preX,f_preY);
        operateDerivative(f_preX,f_preY,*DerivativeMatrixX,*DerivativeMatrixY);
        computeNumericalFlux(f_preX, f_preY);
        operateFlux(*Flux1,*Flux2,*Flux3,*Flux4);
        operateInvereseMass(*MassInverse);
        //q (2) = 0.75*q(0) + 0.25*q(1) + 0.25*∆tR(q (1));
        for(i=0;i<Ney;i++)
            for(j=0;j<Nex;j++)
            {
                cblas_daxpy((N+1)*(N+1),dt,Rate[i][j],1,ConsVariable[i][j],1);
                cblas_dscal((N+1)*(N+1),0.25,ConsVariable[i][j],1);
                cblas_daxpy((N+1)*(N+1),0.75,q0[i][j],1,ConsVariable[i][j],1);
            }

    /*
        computeFlux(f_preX,f_preY);
        operateDerivative(f_preX,f_preY,*DerivativeMatrixX,*DerivativeMatrixY);
        computeNumericalFlux(f_preX, f_preY);
        operateFlux(*Flux1,*Flux2,*Flux3,*Flux4);
        operateInvereseMass(*MassInverse);
        //q (3) = (1/3)*q(0) + (2/3)*q(2) + (2/3)*∆tR(q (2));
        for(i=0;i<Ney;i++)
            for(j=0;j<Nex;j++)
            {
                cblas_daxpy((N+1)*(N+1),dt,Rate[i][j],1,ConsVariable[i][j],1);
                cblas_dscal((N+1)*(N+1),(2.0/3),ConsVariable[i][j],1);
                cblas_daxpy((N+1)*(N+1),(1/3),q0[i][j],1,ConsVariable[i][j],1);
            }
    */
    }

    return ;
}

#endif
