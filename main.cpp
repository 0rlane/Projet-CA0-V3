// DUPORT Orlane & TROUILLARD Thomas
// 21 / 12 / 18

// Projet de CAO :
// Interpolation par élément fini de Powell et Sabin pour la représentation
// d'une surface de classe C1

#include "Point.h"          // classe Point
#include "Triangle.h"       // fonctions d'interpolation
#include "Matrice.h"        // fonctions templates sur matrices
#include "Fonctions_F_G.h"  // fonctions d'interpolations utilisées

#include <iostream>
#include <fstream>
using namespace std;

int main(){

    int NbPts;  // nombre de points dans le domaine D
    int NbTri;  // nombre de triangles dans le domaine D

    ///////////////////////// GENERATION DU MODELE ////////////////////////////////////////

    // ListPoints : tous les points du domaine D (Point NbPts)
    // NT :         les n° de sommets de chaque triangle (int NbTri*3)
    // NTV :        les n° de triangles voisins de chaque triangle (int NbTri*3)
    //                  la valeur -1 est affectée s'il n'y a pas de voisin
    // Omega :      les centre de cercle inscrit de chaque triangle (Point NbTri)
    // NM :         les points Mi associés à chaque triangle (Point NbTri*3)
    // SMT :        sommets des micro-triangles de chaque triangle (Point NbTri*6*3)

    Point*   ListPoints = LecPoints("points.pts",NbPts);
    int**    NT         = lectTriangles("listri.dat",NbTri);
    int**    NTV        = initNTV(NbTri,NT);
    Point*   Omega      = initOmega(ListPoints,NbTri,NT);
    Point**  NM         = initNM(NT,NTV,ListPoints,NbTri,Omega);
    Point*** SMT        = ComputeAllSMT(NbTri, NT, ListPoints, Omega, NM);

    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////// ECRITURE FICHIER RESULTATS //////////////////////////////////

    ofstream fichier("PS.RES"); // fichier de résultats

    results_ListTriangles(fichier, NbTri, NT, Omega); 
    results_ListPoints(fichier, NbTri, NM);

    // output pour chaque fonction étudiée (défini dans Fonctions_F_G.h)
    //      1) fonction exponentielle
    //      2) fonction polynomiale
    for (int NumFonc = 1; NumFonc <= 2; ++NumFonc)
    {
        double** AllCoeff = ComputeAllCoeff(NbTri, NumFonc, NT, ListPoints, NM, Omega);
        // AllCoeff : matrice de tous les coefficients de chaque triangle

        results_ValFonc(fichier, NumFonc, NbPts, ListPoints);
        results_Interpol(fichier, NbTri, NumFonc, ListPoints, NT, Omega, NM, AllCoeff, SMT);
        results_Erreur(fichier, NumFonc, ListPoints, NbPts, NT, Omega, NM, NbTri, AllCoeff, SMT);

        FreeMat(AllCoeff,NbTri);
    }

    fichier.close();

    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////// GENERATION DES SURFACES /////////////////////////////////////

    bool OutputGrille(true);
    // true pour générer les grilles de l'interpolant et des fonctions pour représentation Matlab

    if (OutputGrille)
    {
        ofstream grille1("grille1.pts");
        ofstream fonc1("f1.pts");

        double** AllCoeff;

        AllCoeff = ComputeAllCoeff(NbTri, 1, NT, ListPoints, NM, Omega);
        FichierSurface("grille1.pts","f1.pts",1,ListPoints,NbPts, NT,Omega,NM,NbTri,AllCoeff,SMT);

        AllCoeff = ComputeAllCoeff(NbTri, 2, NT, ListPoints, NM, Omega);
        FichierSurface("grille2.pts","f2.pts",2,ListPoints,NbPts, NT,Omega,NM,NbTri,AllCoeff,SMT);

        FreeMat(AllCoeff,NbTri);
    }

    ///////////////////////////////////////////////////////////////////////////////////////        

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
