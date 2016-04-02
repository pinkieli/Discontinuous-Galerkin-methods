#include <cstdio>
#include <iostream>
#include <string>
#include <cstdlib>
#include "Display.hpp"
using namespace std;
#ifndef Plot_HPP
#define Plot_HPP

void plot(double *X, double *Y, double *Z, unsigned Nx, unsigned Ny, string Name)
{
    int i;
    FILE* pfile;
    pfile=freopen("X.dat","w",stdout);
    display(X,Nx,Ny);
    fclose(stdout);
    pfile=freopen("Y.dat","w",stdout);
    display(Y,Nx,Ny);
    fclose(stdout);
    pfile=freopen("Z.dat","w",stdout);
    display(Z,Nx,Ny);
    fclose(stdout);
    pfile=freopen("temp.m","w",stdout);
    printf("load X.dat;\n");
    printf("load Y.dat;\n");
    printf("load Z.dat;\n");
    printf("surf(X,Y,Z);\n");
    printf("print -djpg %s.jpg;\n",Name.c_str());
    fclose(stdout);

    pfile=freopen ("/dev/tty", "a", stdout);
    i=system("octave temp.m");
    i=system("rm -rf temp* load");

    return ;
}

#endif
