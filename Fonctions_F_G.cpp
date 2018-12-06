#include "Fonctions_F_G.h"
#include <iostream>
#include <math.h>

using namespace std;

double f(Point A){
    // Fonction f(x,y)=exp(x+y)

    double x,y;
    A.getCart(x,y);
    return exp(x+y);
}

double fpx(Point A){
    // Derive en x de la fonction f(x,y)

    double x,y;
    A.getCart(x,y);
    return exp(x+y);
}

double fpy(Point A){
    // Derive en y de la fonction f(x,y)

    double x,y;
    A.getCart(x,y);
    return exp(x+y);
}

Point gradf(Point A){
    // Gradient de la fonction f au point A
    double x,y;
    A.getCart(x,y);
    Point B;
    B.attrib_coord(fpx(A),fpy(A));
    return B;
}

double g(Point A){
    // Fonction g(x,y)

    double x,y;
    A.getCart(x,y);
    return y*y*y-2*x*y*y-5*x*x*y+10*x*y+1;
}

double gpx(Point A){
    // Derive en x de la fonction g(x,y)

    double x,y;
    A.getCart(x,y);
    return -2*y*y-10*x*y+10*y;
}

double gpy(Point A){
    // Derive en y de la fonction g(x,y)

    double x,y;
    A.getCart(x,y);
    return 3*y*y-4*x*y-5*x*x+10*x;
}

Point gradg(Point A){
    // Gradient de la fonction f au point A
    double x,y;
    A.getCart(x,y);
    Point B;
    B.attrib_coord(gpx(A),gpy(A));
    return B;
}

double** initF(int NbPts, Point* ListPoints){
    // Initialisation de la matrice F (NbPts*3) contenant la valeur de f, de fpx et de fpy en chaque point du domaine

    double **F=CreateMat<double>(NbPts,3);

    for (int i=0; i<NbPts; i++){
        F[i][0]=f(ListPoints[i]);  //La premiere colonne contient f(x,y)
        F[i][1]=fpx(ListPoints[i]);  // La seconde colonne contient la derive de f en x
        F[i][2]=fpy(ListPoints[i]);  // La seconde colonne contient la derive de f en y
    }
    return F;
}

double** initG(int NbPts, Point* ListPoints){
    // Initialisation de la matrice F (N*3)

    double **G=CreateMat<double>(NbPts,3);

    for (int i=0; i<NbPts; i++){
        G[i][0]=g(ListPoints[i]);  //La premiere colonne contient g(x,y)
        G[i][1]=gpx(ListPoints[i]);  // La seconde colonne contient la derive de g en x
        G[i][2]=gpy(ListPoints[i]);  // La seconde colonne contient la derive de g en y
    }
    return G;
}

