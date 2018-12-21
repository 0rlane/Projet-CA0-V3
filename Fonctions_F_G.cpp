#include "Fonctions_F_G.h"
#include <iostream>
#include <math.h>

using namespace std;

double f(Point A, int numero_fonction){

    double x,y;
    A.getCart(x,y);

    switch (numero_fonction){
    case 1: // Fonction f(x,y)=exp(x+y)
        return exp(x+y);
    case 2: // Fonction g(x,y)=
        return y*y*y-2*x*y*y-5*x*x*y+10*x*y+1;
    }
}
    
double fpx(Point A, int numero_fonction){

    double x,y;
    A.getCart(x,y);

    switch (numero_fonction){
    case 1: // Fonction f(x,y)=exp(x+y)
        return exp(x+y);
    case 2: // Fonction g(x,y)
        return -2*y*y-10*x*y+10*y;
    }
}

double fpy(Point A, int numero_fonction){

    double x,y;
    A.getCart(x,y);

    switch (numero_fonction){
    case 1: // Fonction f(x,y)=exp(x+y)
        return exp(x+y);
    case 2: // Fonction g(x,y)
        return 3*y*y-4*x*y-5*x*x+10*x;
    }
}

Point gradf(Point A, int numero_fonction){

    double x,y;
    A.getCart(x,y);
    Point B;
    B.attrib_coord(fpx(A,numero_fonction),fpy(A,numero_fonction));

    return B;
}

