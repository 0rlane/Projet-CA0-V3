#include "Point.h"
#include "Triangle.h"
#include "Matrice.h"
#include "Fonctions_F_G.h"
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
    //AfficheMat(NbTri,3,NT);

    /* Creation de la matrice NTV (NbTri*3) contenant le numero de 3 triangles voisins pour chaque triangle
    On affecte la valeur -1 lorsqu'il n'y a pas de triangle voisin */
    int **NTV=initNTV(NbTri,NT);
    //AfficheMat(NbTri,3,NTV);

    // Creation de la matrice Point Omega (NbTri)
    Point *Omega=initOmega(ListPoints,NbTri,NT);
    //afficheCoordPoints(Omega,NbTri);

    // Creation de la matrice Point NM (Nbtri*3)
    Point **NM=initNM(NT,NTV,ListPoints,NbTri,Omega);

    // Calcul des coordonnées barycentriques aux points Mi
    CoordBaryMi(NT, ListPoints, NM, NbTri);
    /*for (int i = 0; i < NbTri; ++i)
    {
        afficheCoordPoints(NM[i], 3);
    }*/

    // Creation du fichier de resultats "PS.RES"
    CreatFileResults("PS.RES",NT,Omega,NM,ListPoints,NbTri,NbPts);

    /*Point A(0.2,1.1);
    int triangle=LocatePointTriangle(A,ListPoints,NT,NbTri);
    cout<<"Le point ";A.affiche();cout<<" est dans le triangle "<<triangle<<endl;
    cout<<"Les sommets de ce triangle sont les points"<<endl;
    cout<<"Point "<<NT[triangle][0]<<" : "; ListPoints[NT[triangle][0]-1].affiche();
    cout<<"Point "<<NT[triangle][1]<<" : "; ListPoints[NT[triangle][1]-1].affiche();
    cout<<"Point "<<NT[triangle][2]<<" : "; ListPoints[NT[triangle][2]-1].affiche();
    int microTriangle=LocatePointMicroTriangle(A,ListPoints,triangle,NT,Omega,NM);
    cout<<"Dans le micro triangle "<<microTriangle<<endl;
    cout<<"Les sommets de ce triangle sont les points"<<endl;
    Point **NMT=MicroTriangle(triangle,NT,ListPoints,Omega,NM);
    cout<<"Point ";NMT[microTriangle][0].affiche();
    cout<<"Point ";NMT[microTriangle][1].affiche();
    cout<<"Point ";NMT[microTriangle][2].affiche();*/


    FreeMat(NT,NbTri); FreeMat(NTV,NbTri);

    delete[] ListPoints;
    delete[] Omega;
    delete[] NM;

    return 0;
}
