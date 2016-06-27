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
#define ABS(a)(a>0.0?a:(-a))
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
    * Lambda    =   The max of absolute value of the eigen value. It is used for compiuting the Rusanov Flux.
    * time      =   The current time step of the Solver.
    */
    double ***eta_RHS,***hu_RHS,***hv_RHS;
    double ***eta_Rate,***hu_Rate,***hv_Rate;
    double ***eta_XFlux,***eta_YFlux,***hu_XFlux,***hu_YFlux,***hv_XFlux,***hv_YFlux,***eta_prev,***hu_prev,***hv_prev;
    double ***eta_XFlux_num,***eta_YFlux_num,***hu_XFlux_num,***hu_YFlux_num,***hv_XFlux_num,***hv_YFlux_num;
    double dx,dy;
    double Lambda;
    double time;

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
    void RK3();
    void updateVelocities();
    void solve();
    void plotSolution(double , double ,string );
    void plotBoundary(double , double , string);
    void writeVTK(string );
};

ShallowWater::ShallowWater(unsigned Nx, unsigned Ny, unsigned n)
{
    /**
    * The constructor for the ShallowWatter class. It takes 3 inputs that defines the solver.
    * The number of elements in the X-direction. [Nex]
    * The number of elements in the Y-direction. [Ney]
    * The order of the polynomials with which each element is to be interpolated.
    */
    Nex =   Nx;
    Ney =   Ny;
    N   =   n;
    time =  0.0;

    /**
    * The memory allocation for all the 3-d/ 2-d array is being done.
    */
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
    /**
    * This function sets the Domain for the solver.
    * It takes 4 inputs :
    * L1: X-starting co-ordinate.
    * L3: Y-starting co-ordinate.
    * L2: X-ending co-ordinate
    * L4: Y-ending co-ordinate.
    */
    unsigned i,j,k1,k2;
    L_start =   L1;
    L_end   =   L2;
    H_start =   L3;
    H_end   =   L4;

    double Xstart,Ystart;
    double Xend,Yend;
    double Nodes[N+1];
    lobattoNodes(Nodes,N+1);

    /**
    * Once the domain is known the `dx` and `dy` can be calculated.
    */
    dx  =   (L_end-L_start)/Nex;
    dy  =   (H_end-H_start)/Ney;

    /**
    * Since, the domain is known, using the information calculating the grid.
    * The variables X[i][j][k] and Y[i][j][k] stores the information about the Domain.
    * Where `i` denotes the serial number of the Element in the Y-direction. Range: [0, Ney)
    * Where `j` denotes the serial number of the Element in the X-direction. Range: [0, Nex)
    * Where `k` denotes the serial number of the grid points in an single Element. Range: [0, (N+1)*(N+1))
    */
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


/**
 * @brief   ShallowWater:: setDepth -- [This function is used to set the Bathymetry `depth` of the Shallow Water. Hence, the membr variable H[i][j][k] is initialized in the function].
 * @param in function<double(double,double)> A [The function takes two values `double x, double y` and returns the Depth. Hence the depth is initialized using a function.]
 */
void ShallowWater::setDepth(function<double(double, double)> A)
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)
    {
        for ( j=0; j<Nex; j++)
        {
            for (k=0;k<((N+1)*(N+1));k++)
                H[i][j][k]  =   A(X[i][j][k],Y[i][j][k]);

            /**
             * Knowing the Bathymetry Height the RHS can be initialized.
             */
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dy,DerivativeMatrixX,(N+1)*(N+1),H[i][j],1,0,H_x[i][j],1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(N+1)*(N+1),(N+1)*(N+1),0.5*dx,DerivativeMatrixY,(N+1)*(N+1),H[i][j],1,0,H_y[i][j],1);
        }
    }

    return ;
}

/**
 * The function `ShallowWater::setInitialConditions(function<double(double,double)> A, function<double(double,double)> B, function<double(double,double)> C)` takes input as 3 functions, each of which separately sets the initial Conditions for the eta(x,y), u(x,y) and v(x,y).
 * @param in function<double(double,double)> A [This function gives the Initial Condition for eta(x,y)]
 * @param in function<double(double,double)> B [This function gives the Initial Condition for u(x,y)]
 * @param in function<double(double,double)> C [This function gives the Initial Condition for v(x,y)]
 */
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

/**
 * The Important parameter of the function i.e. time Step (`dt`) and Number of Time Steps (`NTimeSteps`) is decided in this function.
 * @param a First Inuput is the `dt`
 * @param b Second Input is the `NTimeSteps`
 */
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

/**
 * In each iteration of the RK3, the eta, u, v changes and due to the same point the RHS gets changes, since the `eta` changes.
 */
void ShallowWater::computeRHS()
{
    unsigned i,j,k;
    for(i=0;i<Ney;i++)/*Iterating through the elements along the Y- directiion.*/
        for(j=0;j<Nex;j++)/*Iterating thorugh the elements along the X- direction.*/
            for(k=0;k<((N+1)*(N+1));k++)/*After closing down a single element, iterating through the (N+1)*(N+1) grid points within that element.*/
            {
                eta_RHS[i][j][k]    =   0.0;/*The RHS term for first equation when `eta` is the conservative variable.*/
                hu_RHS[i][j][k]     =   0.0;//-G*eta[i][j][k]*H_x[i][j][k];/*The RHS for the second equation when `hu` is the conservtive variable.*/
                hv_RHS[i][j][k]     =   0.0;//-G*eta[i][j][k]*H_y[i][j][k];/*The RHS term for the third equation when `hv` is the conservative variable.*/
            }
    return ;
}

/**
 * In this function the flux terms for the conservative form of equations is calculated. Each of thd Conservativve variable: eta, eta*u, etau*v has 2 flux terms each. That is the
 * The flux term for the eta is given by \f$\left( \eta u, \eta v\right)\f$.
 * The flux term for the U is given by \f$\left( \eta u^2 + \frac{g*\eta^2}{2}, \eta uv\right)\f$.
 * The flux term for the U is given by \f$\left( \eta uv, \eta v^2 + \frac{g*\eta^2}{2} \right)\f$.
 */
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

/**
 * The numerical flux has to be calculated, which defines the interaction between the adjacent elements. In the following code, Rusanov Flux has been used. The numerical fluz only has to be calculated along the edges. And the same is done.
 */
void ShallowWater::computeNumericalFlux( )
{
    unsigned i,j,k;

    /**
    * The next two big `TRIPLY-NESTED` for loops are for the flux calculations of the internal points.
    */

    //Internal Nodes X-direction.
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

    //Internal Nodes Y-direction
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
    }//Flux- Calculation for the internal nodes end here.


/*
    // Conditions for Periodic Boundary Conditions start here. Hence only the flux at the boundary elements is calculated.
    j=0;
    for(i=0;i<Ney;i++)
    {
        for(k=0;k<=N;k++)
        {
            eta_XFlux_num[i][Nex-1][k*(N+1)+N]  =   eta_XFlux_num[i][j][k*(N+1)]    =   0.5*(eta_XFlux[i][j][k*(N+1)] + eta_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(eta[i][j][k*(N+1)]-eta[i][Nex-1][k*(N+1)+N]) );
            if(ABS(eta_XFlux_num[i][Nex-1][k*(N+1)+N])>1e-8)
                printf("Danda1\n");

            hu_XFlux_num[i][Nex-1][k*(N+1)+N]   =    hu_XFlux_num[i][j][k*(N+1)]    =   0.5*(hu_XFlux[i][j][k*(N+1)] + hu_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(hu[i][j][k*(N+1)]-hu[i][Nex-1][k*(N+1)+N]) );
            if(ABS(hu_XFlux_num[i][Nex-1][k*(N+1)+N]-(hu_XFlux[i][j][k*(N+1)]-(ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)])))*(hu[i][j][k*(N+1)])))>1e-8)
                printf("Danda2\n");

            hv_XFlux_num[i][Nex-1][k*(N+1)+N]   =    hv_XFlux_num[i][j][k*(N+1)]    =   0.5*(hv_XFlux[i][j][k*(N+1)] + hv_XFlux[i][Nex-1][k*(N+1) + N] - (MAX((ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)]))),(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))))*(hv[i][j][k*(N+1)]-hv[i][Nex-1][k*(N+1)+N]) );
            if(ABS(hv_XFlux_num[i][Nex-1][k*(N+1)+N])>1e-8)
                printf("Danda3\n");

        }
    }

    i=0;
    for(j=0;j<Nex;j++)
    {
        for(k=0;k<=N;k++)
        {
            eta_YFlux_num[Ney-1][j][k + N*(N+1)]    =   eta_YFlux_num[i][j][k]  =   0.5*(eta_YFlux[i][j][k] + eta_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(eta[i][j][k] - eta[Ney-1][j][k + N*(N+1)]) );
            if(ABS(eta_YFlux_num[Ney-1][j][k + N*(N+1)])>1e-8)
                printf("Danda4\n");

            hu_YFlux_num[Ney-1][j][k + N*(N+1)]    =    hu_YFlux_num[i][j][k]   =   0.5*(hu_YFlux[i][j][k] + hu_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(hu[i][j][k] - hu[Ney-1][j][k + N*(N+1)]) );
            if(ABS(hu_YFlux_num[Ney-1][j][k + N*(N+1)])>1e-8)
                printf("Danda5\n");

            hv_YFlux_num[Ney-1][j][k + N*(N+1)]    =    hv_YFlux_num[i][j][k]   =   0.5*(hv_YFlux[i][j][k] + hv_YFlux[Ney-1][j][k + N*(N+1)] - (MAX((ABS(v[i][j][k])+sqrt(G*(eta[i][j][k]))),(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))))*(hv[i][j][k] - hv[Ney-1][j][k + N*(N+1)]) );
            if(ABS((hv_YFlux_num[Ney-1][j][k + N*(N+1)])-(hv_YFlux[Ney-1][j][k + N*(N+1)]+(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))*(hv[Ney-1][j][k + N*(N+1)])))>1e-8)
                printf("Final Danda\n");
            if(ABS((hv_YFlux_num[Ney-1][j][k + N*(N+1)])-(hv_YFlux[i][j][k]-(ABS(v[i][j][k])+sqrt(G*(eta[i][j][k])))*(hv[i][j][k])))>1e-8)
                printf("Sorry this is really Final Danda\n");


        }
    } // Periodic Boundary Conditions end here.
*/

   //Writing the Boundary Conditions for the Bath-Tub model.
   //Bath - Tub Model X-flux conditions.
   j=0;
   for(i=0;i<Ney;i++)
   {
       for(k=0;k<=N;k++)
       {
           eta_XFlux_num[i][Nex-1][k*(N+1)+N]   =   0.0 ;
           eta_XFlux_num[i][j][k*(N+1)]         =   0.0 ;

           hu_XFlux_num[i][Nex-1][k*(N+1)+N]    =   hu_XFlux[i][Nex-1][k*(N+1) + N]+(ABS(u[i][Nex-1][k*(N+1) + N])+sqrt(G*(eta[i][Nex-1][k*(N+1) + N])))*(hu[i][Nex-1][k*(N+1)+N]);
           hu_XFlux_num[i][j][k*(N+1)]          =   hu_XFlux[i][j][k*(N+1)]-(ABS(u[i][j][k*(N+1)])+sqrt(G*(eta[i][j][k*(N+1)])))*(hu[i][j][k*(N+1)]);

           hv_XFlux_num[i][Nex-1][k*(N+1)+N]    =   0.0;
           hv_XFlux_num[i][j][k*(N+1)]          =   0.0 ;
       }
   }

   //Bath - Tub Model Y- flux conditions.
   i=0;
   for(j=0;j<Nex;j++)
   {
       for(k=0;k<=N;k++)
       {
           eta_YFlux_num[Ney-1][j][k + N*(N+1)] =   0.0;
           eta_YFlux_num[i][j][k]               =   0.0;

           hu_YFlux_num[Ney-1][j][k + N*(N+1)]  =   0.0;
           hu_YFlux_num[i][j][k]                =   0.0;

           hv_YFlux_num[Ney-1][j][k + N*(N+1)]  =   hv_YFlux[Ney-1][j][k + N*(N+1)]+(ABS(v[Ney-1][j][k + N*(N+1)])+sqrt(G*(eta[Ney-1][j][k + N*(N+1)])))*(hv[Ney-1][j][k + N*(N+1)]);
           hv_YFlux_num[i][j][k]                =   hv_YFlux[i][j][k]-(ABS(v[i][j][k])+sqrt(G*(eta[i][j][k])))*(hv[i][j][k]);
       }
   }

    return ;
}

/**
 * This function operates the Derivative Matrix on the flux terms of the conservative variables, and gets added in the RHS. Here BLAS and LAPACK functions have been used here for performing the linear operations.
 */
void ShallowWater::operateDerivative()
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
/**
 * `ShallowWater::operateFlux` : In this function the Flux Matrices are operated on the Numerical Flux Column vectors, and then added to the RHS.
 */
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

/**
 * `ShallowWater::operateInverseMass` : In this function the Mass Inverse Matrix is operated on the RHS column vector. Hence after operting this, we are only left with the time derivative term on the LHS.
 */

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

/**
 * ShallowWater::copyField : This creates a copy of the conservative variables, as we need to use the earlier time step values at
 */
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
        fprintf(stderr,"Time t=%6.3f\tCourant Number=%6.2f\n",time,Lambda*dt*(1.0/dy+1.0/dx) );
        RK3();
        time += dt;
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

void ShallowWater::plotBoundary(double Y1, double Y2, string outputFilename)
{
    double CGX[Nex*N+1],CGY[Nex*N+1] ;
    zeros(CGX,Nex*N+1);
    zeros(CGY,Nex*N+1);
    unsigned j,k;

    for(j=0;j<Nex;j++)
    {
        k=0;
        CGY[j*N+k]+=0.5*eta[0][j][k];
        CGX[j*N+k]+=0.5*X[0][j][k];
        for(k=1;k<N;k++)
        {
            CGY[j*N+k]+=eta[0][j][k];
            CGX[j*N+k]+=X[0][j][k];
        }

        k=N;
        CGY[j*N+k]+=0.5*eta[0][j][k];
        CGX[j*N+k]+=0.5*X[0][j][k];
    }

    CGY[0]  =   2*CGY[0];
    CGY[Nex*N]=  2*CGY[Nex*N];
    CGX[0]  =   2*CGX[0];
    CGX[Nex*N]=  2*CGX[Nex*N];

    plot(CGX,CGY,Nex*N+1,Y1,Y2,"Shallow Water","Water Height",outputFilename);

    return ;
}

/**
  * Writing the `eta` to VTK files.
  * I am currently going to interpret the field to be of Unstructed_Grid Form of the VTK file formats.
  */
void ShallowWater::writeVTK(string outputFilename)
{
    unsigned i,j,k,k1,k2;
    FILE* pfile;
    pfile   =   freopen(outputFilename.c_str(),"w",stdout);

    // Starting the .vtk files with the cmpulsory field values.
    fprintf(pfile,"# vtk DataFile Version 3.0\nShallow Water\nASCII\nDATASET UNSTRUCTURED_GRID\n");

    // The information of the points.
    fprintf(pfile,"POINTS\t%d\tfloat\n",(N+1)*(N+1)*Nex*Ney);

    // Writing the points and the values.
    for ( i = 0; i < Ney; i++ )
        for ( j = 0; j < Nex; j++ )
            for( k = 0; k < (N+1)*(N+1); k++ )
                fprintf(pfile,"%.3f\t%.3f\t%.3f\n",X[i][j][k],Y[i][j][k],eta[i][j][k]);

    fprintf(pfile,"\n\n");

    // Specifying the information about the CELLS.
    fprintf(pfile,"CELLS\t%u\t%u\n",(N*N*Nex*Ney),5*(N*N*Nex*Ney));

    // Starting writing information about the cells.
    for ( i = 0; i < Ney; i++ )
    {
        for ( j = 0; j < Nex; j++ )
        {
            for( k1 = 0; k1 < N; k1++ )
            {
                for ( k2 = 0; k2 < N; k2++ )
                {
                    k   =   (i*Nex+j)*(N+1)*(N+1) +   k1*(N+1)    +   k2;
                    fprintf(pfile,"%u\t%u\t%u\t%u\t%u\n",4,k,k+1,k+N+2,k+N+1);
                }
            }
        }
    }
    fprintf(pfile,"\n\n");

    // Specifying the information about the CELL TYPES.
    fprintf(pfile,"CELL_TYPES %u\n",(N*N*Nex*Ney));

    // `9` is the CELL TYPE CODE for specifying that it is a quad.
    for ( i = 0; i < (N*N*Nex*Ney); i++)
        fprintf(pfile,"9\n");
    fprintf(pfile,"\n\n");

    // Specifying the information about the POINT_DATA
    fprintf(pfile,"POINT_DATA\t%u\nSCALARS\tEta\tfloat\nLOOKUP_TABLE default\n", (N+1)*(N+1)*Nex*Ney);

    // Writing the value of the POINT_DATA . In this case `eta`.
    for ( i = 0; i < Ney; i++ )
        for ( j = 0; j < Nex; j++ )
            for( k = 0; k < (N+1)*(N+1); k++ )
                fprintf(pfile,"%.3f\n",eta[i][j][k]);

    fclose(pfile);
    return ;
}


#endif
