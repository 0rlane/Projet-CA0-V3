#include "Triangle.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int** lectTriangles(const char* name, int &NbTri){
    /* Recuperation du nombre de triangle (NbTri) dans le fichier name
    et creation de la matrice NT contenant les sommet de chaque traingle (NbTri*3) obtenu dans le fichier name */

    ifstream fichier("listri.dat");
    fichier>>NbTri;

    int** NT=CreateMat<int>(NbTri,3);
    for (int i=0; i<NbTri; i++){
        for (int j=0; j<3; j++){
            fichier>>NT[i][j];
        }
    }
    return NT;
}

int trianglevoisin(int NbTri, int S1, int S2, int k, int **NT){
    /* Recherche du triangle voisin du triangle k passant par les sommets S1 et S2
    Si le traingle n'a pas de voisin on renvoie la valeur -1 */

    int voisin(-1);
    int i(0);
    while (i<NbTri && voisin==-1){
        if (i!=k){
            if ((NT[i][0]==S1 || NT[i][1]==S1 || NT[i][2]==S1) && (NT[i][0]==S2 || NT[i][1]==S2 || NT[i][2]==S2)){
                // Le triangle i est voisin du triangle k
                voisin=i;
            }
        }
    i++;
    }
    return voisin;
}

int** initNTV(int NbTri, int **NT){
    // Initialisation de la matrice NTV (NbTri*3)

    int** NTV=CreateMat<int>(NbTri,3);

    for(int k=0; k<NbTri; k++){
        NTV[k][0]=trianglevoisin(NbTri,NT[k][1],NT[k][2],k,NT);  // triangle voisin des sommets A2A3
        NTV[k][1]=trianglevoisin(NbTri,NT[k][0],NT[k][2],k,NT);  // triangle voisin des sommets A1A3
        NTV[k][2]=trianglevoisin(NbTri,NT[k][0],NT[k][1],k,NT);  // triangle voisin des sommets A1A2
    }
    return NTV;
}

void CalculCoeff(double &a1, double &b1, double &c1, Point A, Point B){
    // Calcul des coefficents a, b, c de la droite (D): ax+by+c=0 passant par les point A et B

    double Xa, Ya, Xb, Yb;
    A.getCoord(Xa,Ya);
    B.getCoord(Xb,Yb);
    b1 = Xb - Xa;
    a1 = Ya - Yb;
    c1 = -Ya*b1 - Xa*a1;
}

Point IntersectionDroites(double a1, double a2, double b1, double b2, double c1, double c2){
    // Calcul les corrdonnées X et Y d'intersection des deux droites de coeff a,b,c

    Point inter;
    double X, Y;

    if (b1!=0. && b2!=0.)
    {
        X = (((c1/b1)-(c2/b2)) / ((a2/b2)-(a1/b1)));
        Y = (-a1/b1)*X-c1/b1;
    }
    else if (b1==0.) // droite 1 est verticale
    {
        X = -c1/a1;
        Y = (-a2/b2)*X-c2/b2;
    }
    else if (b2==0.) // droite 2 est verticale
    {
        X = -c2/a2;
        Y = (-a1/b1)*X-c1/b1;
    }
    else // b1 et b2 ne sont pas sensé etre nuls tous les deux
    cout << "erreur dans IntersectionDroites" << endl;

    inter.attrib_coord(X,Y);
    return inter;
}

double DistEucl(Point A, Point B){
	// Calcule la distance euclidienne entre deux Points

	double Xa, Ya, Xb, Yb;
	A.getCoord(Xa,Ya);
	B.getCoord(Xb,Yb);

	return ( sqrt( (Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb) ) );
}

Point* initOmega(Point* ListPoints, int NbTri, int **NT){
    // Initialisation de la matrice Omega de Nbtri Points

    Point *Omega=new Point[NbTri];

    for (int k=0; k<NbTri; k++){
        double *CentreBary=new double[3];   // Coordonnees barycentrique de omega du triangle k

        // Recuperation des coordonnees des sommets A1, A2 et A3 du triangle k
        Point A1=ListPoints[NT[k][0]-1];
        Point A2=ListPoints[NT[k][1]-1];
        Point A3=ListPoints[NT[k][2]-1];

        double X1,X2,X3,Y1,Y2,Y3;
        A1.getCoord(X1,Y1);
        A2.getCoord(X2,Y2);
        A3.getCoord(X3,Y3);

        // Coordonnees barycentrique de omega du triangle k
        CentreBary[0] = DistEucl(A2,A3);
        CentreBary[1] = DistEucl(A1,A3);
        CentreBary[2] = DistEucl(A1,A2);

        // Coordonnnees cartesiennes de omega du triangle k
        double X, Y;
        X=(CentreBary[0]*X1+CentreBary[1]*X2+CentreBary[2]*X3)/fabs(CentreBary[0]+CentreBary[1]+CentreBary[2]);
        Y=(CentreBary[0]*Y1+CentreBary[1]*Y2+CentreBary[2]*Y3 )/fabs(CentreBary[0]+CentreBary[1]+CentreBary[2]);
        Omega[k].attrib_coord(X,Y);
    }

    return Omega;
}

Point** initNM(int **NT, int **NTV, Point* ListPoints, int nbtri, Point *omega){
    /*Initilisationde la matrice NM Point (nbtriangles*3) qui contient les coordonnees x et y des trois points Mi du triangle k
    NT=matrice (nb_triangles*3) qui contient les 3 sommets pour chaque traingle
    NTV= matrice (nb_triangles*3) qui contient le numero des 3 triangles voisins pour chaque triangle (-1 si il n'y a pas de traingle voisin)
    points=matrice (nb_points*3) qui contient les numero du point et ses coordonnees x et y associees
    nbtri= nombre de triangles dans le domaine
    omega= matrice Point (nb_trinagles) qui contient les coefficients x et y du centre du cercle inscrit de chaque triangle */

    Point** NM=CreateMat<Point>(nbtri,3);

    for (int k=0; k<nbtri; k++){
        double a1, b1, c1;
        double a2, b2, c2;

        for (int i=0; i<3; i++){
            if (NTV[k][i]!=-1){
                // Si il y a un triangle voisin
                CalculCoeff(a1,b1,c1,ListPoints[NT[k][(i+1)%3]-1],ListPoints[NT[k][(i+2)%3]-1]); // Coefficients de la droite (D1): a1*x+b1*y+c1=0 qui passe par les sommets A2 et A3
                CalculCoeff(a2,b2,c2,omega[k],omega[NTV[k][i]]);// Coefficients de la droite (D2): a2*x+b2*y+c2=0 qui passe par les sommets le centre du cercle inscrit du triangle k et par le centre du cercle inscrit du triangle voisin de k par les sommets A2 et A3
                NM[k][i]=IntersectionDroites(a1,a2,b1,b2,c1,c2);

            }else {
                // Si il n'y a pas de triangle voisin
                double X, Y, X1, Y1, X2, Y2;
                ListPoints[NT[k][(i+1)%3]-1].getCoord(X1,Y1);
                ListPoints[NT[k][(i+2)%3]-1].getCoord(X2,Y2);
                X=fabs(X1-X2)/2+min(X1,X2);
                Y=fabs(Y1-Y2)/2+min(Y1,Y2);
                NM[k][i].attrib_coord(X,Y);
            }
        }
    }
    return NM;
}
