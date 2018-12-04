#include "Point.h"
#include <iostream>
#include <fstream>
using namespace std;

/////////////////////////////////////////////// FONCTION POINT /////////////////////////////////////////////////

Point::Point( double xi, double yi):X(xi),Y(yi){}

Point::~Point() {}

void Point::affiche()const{
    /* Affichage des coordonnee du Point */

	cout << X << " " << Y << endl;
}

void Point::attrib_coord(double xi, double yi){
     /* Attribution des coordonnees du Point */

	X = xi;
	Y = yi;
}

void Point::getCart(double &xi,double &yi){
    /* Renvoie la valeur des coordonnees du point */

	xi = X;
	yi = Y;
}

//////////////////////////////////////////////// FONCTIONS NON-MEMBRES ////////////////////////////////////////////

Point* LecPoints(const char* name, int &N){
    /* Lecture du fichier name pour extraire le nombre de points dans le domaine D (N)
    et pour extraire les coordonnees de chaque point mis dans un vecteur de N Points */

	double tmp, A, B;

	ifstream fichier(name);
	fichier >> N;

	Point* listPoints = new Point[N];
	for (int i = 0; i < N; ++i)
	{
		fichier >> tmp >> A >> B;
		listPoints[i].attrib_coord(A,B);
	}

	fichier.close();
	return listPoints;
}

void afficheCoordPoints(Point* ListPoints, int N){
    // Affichage des traingles et de leurs coordonnees cartesiennes

    for(int k=0; k<N; k++){
        cout<<"Triangle "<<k<<": ";
        ListPoints[k].affiche();
    }
}

void CreateFilePoint(const char* name, int nrow, Point *A){
    // Creation d'un fichier qui comprend les coordonnnees des points du vecteur A

    double X,Y;
    ofstream fichier(name);
    for (int i=0; i<nrow; i++){
        A[i].getCoord(X,Y);
        fichier<<X<<" "<<Y<<endl;
    }
    fichier.close();
}

void CreateFileMatPoint(const char* name, int nrow, int ncol, Point **A){
    // Creation d'un fichier qui comprend les coordonnnees des points du vecteur A

    double X,Y;
    ofstream fichier(name);
    for (int i=0; i<nrow; i++){
        for (int j=0; j<ncol; j++){
            A[i][j].getCoord(X,Y);
            fichier<<X<<" "<<Y<<" ";
        }
        fichier<<endl;
    }
    fichier.close();
}

