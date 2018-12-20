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

    ///////////////////////// GENERATION DU MODELE ////////////////////////////////////////

    // Creation d'un vecteur de nbPts Point contenant tous les points du domaine D
    Point* ListPoints=LecPoints("points.pts",NbPts);

    // Creation de la matrice NM (NbTri*3) contenant les 3 sommets de chaque triangle
    int **NT=lectTriangles("listri.dat",NbTri);

    /* Creation de la matrice NTV (NbTri*3) contenant le numero de 3 triangles voisins pour chaque triangle
    On affecte la valeur -1 lorsqu'il n'y a pas de triangle voisin */
    int **NTV=initNTV(NbTri,NT);

    // Creation de la matrice Point Omega (NbTri)
    Point *Omega=initOmega(ListPoints,NbTri,NT);

    // Creation de la matrice Point NM (Nbtri*3)
    Point **NM=initNM(NT,NTV,ListPoints,NbTri,Omega);

    // Calcul des coordonnées barycentriques aux points Mi
    CoordBaryMi(NT, ListPoints, NM, NbTri);

    // Calcule des sommets des MicroTriangle pour chaque triangle
    Point*** SMT = ComputeAllSMT(NbTri, NT, ListPoints, Omega, NM);

    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////// ECRITURE FICHIER RESULTATS //////////////////////////////////

    // output
    ofstream fichier("PS.RES");

    results_ListTriangles(fichier, NbTri, NT, Omega);
    results_ListPoints(fichier, NbTri, NM);

    // output pour chaque fonction étudiée (défini dans Fonctions_F_G.h)
    //      1) fonction exponentielle
    //      2) fonction polynomiale
    for (int NumFonc = 1; NumFonc <= 2; ++NumFonc)
    {
        // Calcule la matrice de tous les coefficients de chaque triangle
        double** AllCoeff = ComputeAllCoeff(NbTri, NumFonc, NT, ListPoints, NM, Omega);

        results_ValFonc(fichier, NumFonc, NbPts, ListPoints);
        results_Interpol(fichier, NbTri, NumFonc, ListPoints, NT, Omega, NM, AllCoeff, SMT);
        results_Erreur(fichier, NumFonc, ListPoints, NbPts, NT, Omega, NM, NbTri, AllCoeff, SMT);

        FreeMat(AllCoeff,NbTri);
    }

    fichier.close();

    ///////////////////////////////////////////////////////////////////////////////////////

    
    double **Grille;
    double** AllCoeff = ComputeAllCoeff(NbTri, 1, NT, ListPoints, NM, Omega);
    Grille=InterpolantDomaine(ListPoints,NbPts, NT,Omega,NM,NbTri,AllCoeff,SMT);
    ofstream grille("grille1.pts");
    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){
            grille<<Grille[i][j]<<" ";
        }
        grille<<endl;
    }
    grille.close();
    FreeMat(AllCoeff,NbTri);
    

    FreeMat(NT,NbTri); FreeMat(NTV,NbTri);

    for (int i = 0; i < NbTri; ++i)
    {
        FreeMat(SMT[i],6);
    }

    delete[] ListPoints;
    delete[] Omega;
    delete[] NM;

    return 0;
}
