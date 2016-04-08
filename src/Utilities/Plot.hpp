#include <cstdio>
#include <iostream>
#include <string>
#include <cstdlib>
#include "Display.hpp"
using namespace std;
#ifndef Plot_HPP
#define Plot_HPP

void plot(double *X, double *Y, double *Z, unsigned Nx, unsigned Ny,double X1, double X2,double Y1, double Y2,double Z1, double Z2, string Name)
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
    printf("h=figure(1);\n");
    printf("surf(X,Y,Z);\n");
    printf("axis([%.2f %.2f %.2f %.2f %.2f %.2f]);\n",X1,X2,Y1,Y2,Z1,Z2);
    printf("H = 8; W = 8;\nset(h,\'PaperUnits\',\'inches\');\nset(h,\'PaperOrientation\',\'portrait\');\nset(h,\'PaperSize\',[H,W]);\n");
    printf("print(h,\'-djpg\',\'-color\',\'%s.jpg\');\n",Name.c_str());
    fclose(stdout);

    pfile=freopen ("/dev/tty", "a", stdout);
    i=system("octave temp.m");
    i=system("rm -rf temp* load");

    return ;
}

#endif
