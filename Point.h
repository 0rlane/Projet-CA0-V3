#ifndef POINT_H
#define POINT_H


class Point{
private:
    // Donnees
    double X,Y; // Coordonnees cartesiennes
    double w1,w2,w3; //coordonnees barycentriques

public:
    // constructeur & destructeur
	Point() { X = Y = -1.; }
	Point( double xi, double yi);
	~Point();

	// fonctions
	void affiche(void)const;
	void attrib_coord(double x, double y);
	void attrib_bary(double c1, double c2, double c3);
    void getCart(double &x,double &y);
    void getBary(double &c1, double &c2, double &c3);

};

// Fonctions externes mais qui sont en rapport avec la classe Point
Point* LecPoints(const char* name, int &N);
void afficheCoordPoints(Point* ListPoints, int N);
void CreateFilePoint(const char* name, int nrow, Point *A);
void CreateFileMatPoint(const char* name, int nrow, int ncol, Point **A);
double multPointsCart(Point A, Point B);
Point calcVect(Point A, Point B);

#endif // POINT_H

