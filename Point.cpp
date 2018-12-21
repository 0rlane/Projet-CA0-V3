#include "Point.h"
#include <iostream>
#include <fstream>
using namespace std;

/////////////////////////////////////////////// FONCTION POINT /////////////////////////////////////////////////

Point::Point( double xi, double yi):X(xi),Y(yi){}

Point::~Point() {}

void Point::affiche()const{
    /* Affichage des coordonnee du Point */

	cout << "cartesiennes : "<< X << " " << Y << " barycentriques : " << w1 << " " << w2 << " " << w3 << endl;
}

void Point::attrib_coord(double xi, double yi){
     /* Attribution des coordonnees cartesiennes du Point */

	X = xi;
	Y = yi;
}

void Point::attrib_bary(double c1, double c2, double c3){
     /* Attribution des coordonnees barycentriques du Point */
    w1 = c1;
    w2 = c2;
    w3 = c3;
}

void Point::getCart(double &xi,double &yi){
    /* Renvoie la valeur des coordonnees cartesiennes du Point */

	xi = X;
	yi = Y;
}

void Point::getBary(double &c1, double &c2, double &c3){
    /* Renvoie la valeur des coordonnees barycentriques du Point */
    c1 = w1;
    c2 = w2;
    c3 = w3;
 }

//////////////////////////////////////////////// FONCTIONS NON-MEMBRES ////////////////////////////////////////////

Point* LecPoints(const char* name, int &N){
    /* Lecture du fichier "name" et extraction du nombre de points et de leurs coordonnées 
    Renvoie le vecteur de Point (Point NbPts)*/

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

double multPointsCart(Point A, Point B){
    // Multiplication des coordonnees cartesiennes des point A et B

    double xA,yA,xB,yB;
    A.getCart(xA,yA);
    B.getCart(xB,yB);

    double C;
    C=xA*xB+yA*yB;
    return C;
}

Point calcVect(Point A, Point B){
    // Calcule des composantes du vecteur AB

    double xA,yA,xB,yB;
    A.getCart(xA,yA);
    B.getCart(xB,yB);
    Point C(xB-xA,yB-yA);
    return C;
}
