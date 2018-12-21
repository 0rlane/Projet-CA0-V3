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

// Fonctions externes liées à la classe Point
Point* LecPoints(const char* name, int &N);
double multPointsCart(Point A, Point B);
Point calcVect(Point A, Point B);

#endif // POINT_H

