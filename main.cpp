#include "Point.h"
#include "Triangle.h"
#include "Matrice.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
    int NbPts;  // nombre de points dans le domaine D
    int NbTri;  // nombre de triangles dans le domaine D

    // Creation d'un vecteur de nbPts Point contenant tous les points du domaine D
    Point* ListPoints=LecPoints("points.pts",NbPts);
    //afficheCoordPoints(ListPoints,NbPts);

    // Creation de la matrice NM (NbTri*3) contenant les 3 sommets de chaque triangle
    int **NT=lectTriangles("listri.dat",NbTri);
    //AfficheMat(NbTri,3,NM);

    /* Creation de la matrice NTV (NbTri*3) contenant le numero de 3 triangles voisins pour chaque triangle
    On affecte la valeur -1 lorsqu'il n'y a pas de triangle voisin */
    int **NTV=initNTV(NbTri,NT);
    //AfficheMat(NbTri,3,NTV);

    // Creation de la matrice Point Omega (NbTri)
    Point *Omega=initOmega(ListPoints,NbTri,NT);
    //afficheCoordPoints(Omega,NbTri);

    // Creation de la matrice Point NM (Nbtri*3)
    Point **NM=initNM(NT,NTV,ListPoints,NbTri,Omega);




    FreeMat(NT,NbTri); FreeMat(NTV,NbTri);

    return 0;
}

