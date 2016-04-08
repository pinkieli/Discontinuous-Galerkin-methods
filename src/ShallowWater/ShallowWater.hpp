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
#define G  9.81

using namespace std;

#ifndef ShallowWater_HPP
#define ShallowWater_HPP

class ShallowWater
{
private:

    /**
    * Defining the solver properties.
    * Nex   =   The  number of elements used in the x-direction.
    * Ney   =   The number of elements used in the y-direction.
    * N     =   The order of polynomial of the basis function used while approaching the approximate solution.
    */
    unsigned Nex,Ney,N;

    /**
    * Currently the domain is assumed to be a rectangle.
    * L_start   =   The x-starting co-ordinate.
    * H_start   =   The y-starting co-ordinate.
    * L_end     =   The x-ending co-ordinate.
    * H_end     =   The y-ending co-ordinate.
    */
    double L_start,L_end,H_start,H_end;

    /**
    * Nomenclature:
    * eta   =   The total fluid height from the bottom.
    * H     =   The total height of the ground/ mean level of the fluid.
    * H_x   =   Partial Derivative of `H` wrt X.
    * H_y   =   Partial Derivative of `H` wrt X.
    * u     =   The veolcity in the x-direction.
    * v     =   The velocity in the y-direction.
    * hu    =   The second term on which conservation is applied.
    * hv    =   The third term on which the conservation is applied.
    */
    double ***eta,***H,***H_x,***H_y,***u,***v,***hu,***hv;

    /**
    * X =   The x-co-ordinate of the grid points.
    * Y =   The y- co-ordinate of the grid points.
    */
    double ***X, ***Y;

    /**
    * Solver Properties.
    * dt            =    The time step.
    * NTimeSteps    =   The number of timesteps needed for calculating the solution.
    */
    double dt;
    unsigned NTimeSteps;

    /**
    * In- class created variables the variables ahead have no physical significance, just are defined as function properties so as to faciliate the transmission of data from one member function to the another.
    *  Lambda   =   The max of absolute value of the eigen value. It is used for compiuting the Rusanov Flux.
    */
    double ***eta_RHS,***hu_RHS,***hv_RHS;
    double ***eta_Rate,***hu_Rate,***hv_Rate;
    double ***eta_XFlux,***eta_YFlux,***hu_XFlux,***hu_YFlux,***hv_XFlux,***hv_YFlux,***eta_prev,***hu_prev,***hv_prev;
    double ***eta_XFlux_num,***eta_YFlux_num,***hu_XFlux_num,***hu_YFlux_num,***hv_XFlux_num,***hv_YFlux_num;
    double dx,dy;
    double Lambda;

    /**
    *   The matrices which will be used for solution methodology purposes.
    */
    double *MassMatrix, *MassInverse, *DerivativeMatrixX, *DerivativeMatrixY,*Flux1,*Flux2,*Flux3,*Flux4;

public:

    ShallowWater(unsigned , unsigned, unsigned );

    void setDomain(double , double , double , double );
    void setDepth(function<double(double,double)> );
    void setInitialConditions(function<double(double,double)> , function<double(double,double)>, function<double(double,double)> );
    void setSolver(double , unsigned );
    void computeLambda();
    void computeFlux();
    void computeRHS();
    void operateDerivative();
    void computeNumericalFlux();
    void operateFlux();
    void operateInverseMass();
    void copyField();
    void plotSolution(double , double ,string );
    void RK3();
    void updateVelocities();
    void solve();
};

ShallowWater::ShallowWater(unsigned Nx, unsigned Ny, unsigned n)
{
    Nex =   Nx;
    Ney =   Ny;
    N   =   n;

    eta             =   create3D(Ney,Nex,(N+1)*(N+1));
    H               =   create3D(Ney,Nex,(N+1)*(N+1));
    H_x             =   create3D(Ney,Nex,(N+1)*(N+1));
    H_y             =   create3D(Ney,Nex,(N+1)*(N+1));
    u               =   create3D(Ney,Nex,(N+1)*(N+1));
    v               =   create3D(Ney,Nex,(N+1)*(N+1));
    hu              =   create3D(Ney,Nex,(N+1)*(N+1));
    hv              =   create3D(Ney,Nex,(N+1)*(N+1));
    X               =   create3D(Ney,Nex,(N+1)*(N+1));
    Y               =   create3D(Ney,Nex,(N+1)*(N+1));

    eta_RHS         =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_RHS          =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_RHS          =   create3D(Ney,Nex,(N+1)*(N+1));

    eta_Rate        =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_Rate         =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_Rate         =   create3D(Ney,Nex,(N+1)*(N+1));

    eta_XFlux       =   create3D(Ney,Nex,(N+1)*(N+1));
    eta_YFlux       =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_XFlux        =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_YFlux        =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_XFlux        =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_YFlux        =   create3D(Ney,Nex,(N+1)*(N+1));
    eta_XFlux_num   =   create3D(Ney,Nex,(N+1)*(N+1));
    eta_YFlux_num   =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_XFlux_num    =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_YFlux_num    =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_XFlux_num    =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_YFlux_num    =   create3D(Ney,Nex,(N+1)*(N+1));

    eta_prev           =   create3D(Ney,Nex,(N+1)*(N+1));
    hu_prev            =   create3D(Ney,Nex,(N+1)*(N+1));
    hv_prev            =   create3D(Ney,Nex,(N+1)*(N+1));

    MassMatrix          =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    MassInverse         =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    DerivativeMatrixX   =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    DerivativeMatrixY   =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    Flux1               =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    Flux2               =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    Flux3               =   new double[((N+1)*(N+1))*((N+1)*(N+1))];
    Flux4               =   new double[((N+1)*(N+1))*((N+1)*(N+1))];

    twoDMassMatrix(MassMatrix,N);
    inverse(MassMatrix,MassInverse,(N+1)*(N+1));

    twoDDerivativeMatrixX(DerivativeMatrixX,N);
    twoDDerivativeMatrixY(DerivativeMatrixY,N);

    twoDFluxMatrix1(Flux1,N);
    twoDFluxMatrix2(Flux2,N);
    twoDFluxMatrix3(Flux3,N);
    twoDFluxMatrix4(Flux4,N);
}

void ShallowWater::setDomain(double L1, double L2, double L3, double L4)
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


void ShallowWater::setDepth(function<double(double, double)> A)
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for ( j=0; j<Nex; j++)
        {
            for (k=0;k<((N+1)*(N+1));k++)
                H[i][j][k]  =   A(X[i][j][k],Y[i][j][k]);

            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),H[i][j],1,0,H_x[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),H[i][j],1,0,H_y[i][j],1);
        }
    }

    return ;
}


void ShallowWater::setInitialConditions(function<double(double,double)> A, function<double(double,double)> B, function<double(double,double)> C)
{
    unsigned i,j,k;

    for ( i=0;i<Ney;i++)
        for ( j=0; j<Nex; j++)
            for (k=0;k<((N+1)*(N+1));k++)
            {
                eta[i][j][k]    =   A(X[i][j][k],Y[i][j][k])  +   H[i][j][k];
                u[i][j][k]      =   B(X[i][j][k],Y[i][j][k]);
                v[i][j][k]      =   C(X[i][j][k],Y[i][j][k]);

                hu[i][j][k]     =   eta[i][j][k]*u[i][j][k];
                hv[i][j][k]     =   eta[i][j][k]*v[i][j][k];
            }

    return ;
}

void ShallowWater::setSolver(double a, unsigned b)
{
    dt          =   a;
    NTimeSteps  =   b;
}

void ShallowWater::computeLambda()
{

    Lambda = 0;
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                Lambda  =   MAX(Lambda,ABS(v[i][j][k]));
                Lambda  =   MAX(Lambda,ABS(u[i][j][k*(N+1)]));
            }
        }
    }

    j   =   Nex-1;

    for(i=0;i<Ney;i++)
        for(k=0;k<=N;k++)
            Lambda  =   MAX(Lambda,ABS(u[i][j][k*(N+1)+N]));

    i   =   Ney-1;

    for(j=0;j<Nex;j++)
        for(k=0;k<=N;k++)
            Lambda  =   MAX(Lambda,ABS(v[i][j][k + N*(N+1) ]));

    return ;
}

void ShallowWater::computeRHS()
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
            for(k=0;k<((N+1)*(N+1));k++)
            {
                eta_RHS[i][j][k]    =   0.0;
                hu_RHS[i][j][k]     =   -G*eta[i][j][k]*H_x[i][j][k];
                hv_RHS[i][j][k]     =   -G*eta[i][j][k]*H_y[i][j][k];
                //hu_RHS[i][j][k]     =   0.0;
                //hv_RHS[i][j][k]     =   0.0;
            }


    return ;
}


void ShallowWater::computeFlux()
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<((N+1)*(N+1));k++)
            {
                eta_XFlux[i][j][k]  = eta[i][j][k]*u[i][j][k];
                eta_YFlux[i][j][k]  = eta[i][j][k]*v[i][j][k];

                hu_XFlux[i][j][k]  = eta[i][j][k]*u[i][j][k]*u[i][j][k] +   0.5*G*eta[i][j][k]*eta[i][j][k];
                hu_YFlux[i][j][k]  = eta[i][j][k]*u[i][j][k]*v[i][j][k];

                hv_XFlux[i][j][k]  = eta[i][j][k]*u[i][j][k]*v[i][j][k];
                hv_YFlux[i][j][k]  = eta[i][j][k]*v[i][j][k]*v[i][j][k] +   0.5*G*eta[i][j][k]*eta[i][j][k];
            }
        }
    }
    return ;
}
void ShallowWater::computeNumericalFlux( )
{
    unsigned i,j,k;

    for(i=0;i<Ney;i++)
    {
        for(j=1;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                eta_XFlux_num[i][j-1][k*(N+1)+N]    =   eta_XFlux_num[i][j][k*(N+1)]    =   0.5*(eta_XFlux[i][j][k*(N+1)] + eta_XFlux[i][j-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][j-1][k*(N+1)+N])+sqrt(G*(eta[i][j-1][k*(N+1)+N])))))*(eta[i][j][k*(N+1)]-eta[i][j-1][k*(N+1)+N]) );

                hu_XFlux_num[i][j-1][k*(N+1)+N]     =   hu_XFlux_num[i][j][k*(N+1)]    =   0.5*(hu_XFlux[i][j][k*(N+1)] + hu_XFlux[i][j-1][k*(N+1) + N] - MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][j-1][k*(N+1)+N])+sqrt(G*(eta[i][j-1][k*(N+1)+N]))))*(hu[i][j][k*(N+1)]-hu[i][j-1][k*(N+1)+N]) );

                hv_XFlux_num[i][j-1][k*(N+1)+N]     =   hv_XFlux_num[i][j][k*(N+1)]    =   0.5*(hv_XFlux[i][j][k*(N+1)] + hv_XFlux[i][j-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][j-1][k*(N+1)+N])+sqrt(G*(eta[i][j-1][k*(N+1)+N])))))*(hv[i][j][k*(N+1)]-hv[i][j-1][k*(N+1)+N]) );
            }
        }
    }

    j=0;
    for(i=0;i<Ney;i++)
    {
        for(k=0;k<=N;k++)
        {
            eta_XFlux_num[i][Nex-1][k*(N+1)+N]  =   eta_XFlux_num[i][j][k*(N+1)]    =   0.5*(eta_XFlux[i][j][k*(N+1)] + eta_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(eta[i][j][k*(N+1)]-eta[i][Nex-1][k*(N+1)+N]) );

            hu_XFlux_num[i][Nex-1][k*(N+1)+N]   =    hu_XFlux_num[i][j][k*(N+1)]    =   0.5*(hu_XFlux[i][j][k*(N+1)] + hu_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(hu[i][j][k*(N+1)]-hu[i][Nex-1][k*(N+1)+N]) );

            hv_XFlux_num[i][Nex-1][k*(N+1)+N]   =    hv_XFlux_num[i][j][k*(N+1)]    =   0.5*(hv_XFlux[i][j][k*(N+1)] + hv_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(hv[i][j][k*(N+1)]-hv[i][Nex-1][k*(N+1)+N]) );

        }
    }

    for(i=1;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            for(k=0;k<=N;k++)
            {
                eta_YFlux_num[i-1][j][k + N*(N+1)]  =   eta_YFlux_num[i][j][k]  =   0.5*(eta_YFlux[i][j][k] + eta_YFlux[i-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[i-1][j][k + N*(N+1)])+sqrt(G*(eta[i-1][j][k + N*(N+1)])))))*(eta[i][j][k] - eta[i-1][j][k + N*(N+1)]) );

                hu_YFlux_num[i-1][j][k + N*(N+1)]  =    hu_YFlux_num[i][j][k]   =   0.5*(hu_YFlux[i][j][k] + hu_YFlux[i-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[i-1][j][k + N*(N+1)])+sqrt(G*(eta[i-1][j][k + N*(N+1)])))))*(hu[i][j][k] - hu[i-1][j][k + N*(N+1)]) );

                hv_YFlux_num[i-1][j][k + N*(N+1)]  =    hv_YFlux_num[i][j][k]   =   0.5*(hv_YFlux[i][j][k] + hv_YFlux[i-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[i-1][j][k + N*(N+1)])+sqrt(G*(eta[i-1][j][k + N*(N+1)])))))*(hv[i][j][k] - hv[i-1][j][k + N*(N+1)]) );
            }
        }
    }

    i=0;

    for(j=0;j<Nex;j++)
    {
        for(k=0;k<=N;k++)
        {
            eta_YFlux_num[Ney-1][j][k + N*(N+1)]    =   eta_YFlux_num[i][j][k]  =   0.5*(eta_YFlux[i][j][k] + eta_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(eta[i][j][k] - eta[Ney-1][j][k + N*(N+1)]) );

            hu_YFlux_num[Ney-1][j][k + N*(N+1)]    =    hu_YFlux_num[i][j][k]   =   0.5*(hu_YFlux[i][j][k] + hu_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(hu[i][j][k] - hu[Ney-1][j][k + N*(N+1)]) );

            hv_YFlux_num[Ney-1][j][k + N*(N+1)]    =    hv_YFlux_num[i][j][k]   =   0.5*(hv_YFlux[i][j][k] + hv_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(hv[i][j][k] - hv[Ney-1][j][k + N*(N+1)]) );

        }
    }

    return ;
}
void ShallowWater::operateDerivative( )
{
    unsigned i,j;

    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),eta_XFlux[i][j],1,1,eta_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),eta_YFlux[i][j],1,1,eta_RHS[i][j],1);

            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),hu_XFlux[i][j],1,1,hu_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),hu_YFlux[i][j],1,1,hu_RHS[i][j],1);

            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),hv_XFlux[i][j],1,1,hv_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),hv_YFlux[i][j],1,1,hv_RHS[i][j],1);
        }
    }
    return ;
}

void ShallowWater::operateFlux()
{
    unsigned i,j;

    for(i=0;i<Ney;i++)
    {
        for(j=0;j<Nex;j++)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dy,Flux2,(N+1)*(N+1),eta_XFlux_num[i][j],1,1,eta_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dy,Flux4,(N+1)*(N+1),eta_XFlux_num[i][j],1,1,eta_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dx,Flux3,(N+1)*(N+1),eta_YFlux_num[i][j],1,1,eta_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dx,Flux1,(N+1)*(N+1),eta_YFlux_num[i][j],1,1,eta_RHS[i][j],1);

            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dy,Flux2,(N+1)*(N+1),hu_XFlux_num[i][j],1,1,hu_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dy,Flux4,(N+1)*(N+1),hu_XFlux_num[i][j],1,1,hu_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dx,Flux3,(N+1)*(N+1),hu_YFlux_num[i][j],1,1,hu_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dx,Flux1,(N+1)*(N+1),hu_YFlux_num[i][j],1,1,hu_RHS[i][j],1);

            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dy,Flux2,(N+1)*(N+1),hv_XFlux_num[i][j],1,1,hv_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dy,Flux4,(N+1)*(N+1),hv_XFlux_num[i][j],1,1,hv_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),-0.5*dx,Flux3,(N+1)*(N+1),hv_YFlux_num[i][j],1,1,hv_RHS[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1), 0.5*dx,Flux1,(N+1)*(N+1),hv_YFlux_num[i][j],1,1,hv_RHS[i][j],1);
        }
    }

    return ;
}

void ShallowWater::operateInverseMass()
{
    unsigned i,j;
    double alpha    =   4.0/(dx*dy);
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),alpha,MassInverse,(N+1)*(N+1),eta_RHS[i][j],1,0,eta_Rate[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),alpha,MassInverse,(N+1)*(N+1),hu_RHS[i][j],1,0,hu_Rate[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),alpha,MassInverse,(N+1)*(N+1),hv_RHS[i][j],1,0,hv_Rate[i][j],1);
        }


    return ;
}

void ShallowWater::copyField()
{
    unsigned i,j;
    for( i = 0; i< Ney; i++)
        for(j=0;j<Nex;j++)
        {
            memcpy(eta_prev[i][j],eta[i][j],(N+1)*(N+1)*sizeof(double));
            memcpy(hu_prev[i][j],hu[i][j],(N+1)*(N+1)*sizeof(double));
            memcpy(hv_prev[i][j],hv[i][j],(N+1)*(N+1)*sizeof(double));
        }

    return ;
}

void ShallowWater::updateVelocities()
{
    unsigned i,j,k;

    for ( i=0;i<Ney;i++)
        for ( j=0; j<Nex; j++)
            for (k=0;k<((N+1)*(N+1));k++)
            {
                u[i][j][k]      =   hu[i][j][k]/eta[i][j][k];
                v[i][j][k]      =   hv[i][j][k]/eta[i][j][k];
            }
    return ;
}

void ShallowWater::RK3()
{
    unsigned i,j;
    computeLambda();
    copyField();

    computeRHS();
    computeFlux();
    computeNumericalFlux();
    operateDerivative();
    operateFlux();
    operateInverseMass();
    //q (1) = q (0) + ∆tR(q (0));
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
        {
            cblas_daxpy((N+1)*(N+1),dt,eta_Rate[i][j],1,eta[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hu_Rate[i][j],1,hu[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hv_Rate[i][j],1,hv[i][j],1);
        }
    updateVelocities();

    computeRHS();
    computeFlux();
    computeNumericalFlux();
    operateDerivative();
    operateFlux();
    operateInverseMass();
    //q (2) = 0.5*q(0) + 0.5*q(1) + 0.5*∆tR(q (1));
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
        {
            cblas_daxpy((N+1)*(N+1),dt,eta_Rate[i][j],1,eta[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hu_Rate[i][j],1,hu[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hv_Rate[i][j],1,hv[i][j],1);
            cblas_dscal((N+1)*(N+1),0.25,eta[i][j],1);
            cblas_dscal((N+1)*(N+1),0.25,hu[i][j],1);
            cblas_dscal((N+1)*(N+1),0.25,hv[i][j],1);
            cblas_daxpy((N+1)*(N+1),0.75,eta_prev[i][j],1,eta[i][j],1);
            cblas_daxpy((N+1)*(N+1),0.75,hu_prev[i][j],1,hu[i][j],1);
            cblas_daxpy((N+1)*(N+1),0.75,hv_prev[i][j],1,hv[i][j],1);
        }
    updateVelocities();

    computeRHS();
    computeFlux();
    computeNumericalFlux();
    operateDerivative();
    operateFlux();
    operateInverseMass();
    //q (2) = 0.5*q(0) + 0.5*q(1) + 0.5*∆tR(q (1));
    for(i=0;i<Ney;i++)
        for(j=0;j<Nex;j++)
        {
            cblas_daxpy((N+1)*(N+1),dt,eta_Rate[i][j],1,eta[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hu_Rate[i][j],1,hu[i][j],1);
            cblas_daxpy((N+1)*(N+1),dt,hv_Rate[i][j],1,hv[i][j],1);
            cblas_dscal((N+1)*(N+1),(2.0/3.0),eta[i][j],1);
            cblas_dscal((N+1)*(N+1),(2.0/3.0),hu[i][j],1);
            cblas_dscal((N+1)*(N+1),(2.0/3.0),hv[i][j],1);
            cblas_daxpy((N+1)*(N+1),(1.0/3.0),eta_prev[i][j],1,eta[i][j],1);
            cblas_daxpy((N+1)*(N+1),(1.0/3.0),hu_prev[i][j],1,hu[i][j],1);
            cblas_daxpy((N+1)*(N+1),(1.0/3.0),hv_prev[i][j],1,hv[i][j],1);
        }
    updateVelocities();


    return ;
}

void ShallowWater::solve()
{
    unsigned t;

    for(t=0;t<NTimeSteps;t++)
    {
        printf("Time t=%6.3f\tCourant Number=%6.2f\n",(t+1.0)*dt,Lambda*dt*(1.0/dy+1.0/dx) );
        RK3();
    }

    return ;
}

void ShallowWater::plotSolution(double Z1, double Z2, string s)
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
                CG[i*N+k1][j*N+k2]+=0.5*eta[i][j][k1*(N+1)+k2];

                for(k2=1;k2<N;k2++)
                    CG[i*N+k1][j*N+k2]=eta[i][j][k1*(N+1)+k2];

                k2=N;
                CG[i*N+k1][j*N+k2]+=0.5*eta[i][j][k1*(N+1)+k2];
            }

            k1=0;
            k2=0;
            CG[i*N+k1][j*N+k2]+=0.25*eta[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CG[i*N+k1][j*N+k2]+=0.5*eta[i][j][k1*(N+1)+k2];
            k2=N;
            CG[i*N+k1][j*N+k2]+=0.25*eta[i][j][k1*(N+1)+k2];

            k1=N;
            k2=0;
            CG[i*N+k1][j*N+k2]+=0.25*eta[i][j][k1*(N+1)+k2];
            for(k2=1;k2<N;k2++)
                CG[i*N+k1][j*N+k2]+=0.5*eta[i][j][k1*(N+1)+k2];
            k2=N;
            CG[i*N+k1][j*N+k2]+=0.25*eta[i][j][k1*(N+1)+k2];
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

    plot(*CGX,*CGY,*CG,Ney*N+1,Nex*N+1,L_start,L_end,H_start,H_end,Z1,Z2,s);

    return ;
}


#endif
