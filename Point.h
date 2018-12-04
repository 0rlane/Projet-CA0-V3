#ifndef POINT_H
#define POINT_H


class Point{
private:
    // Donnees
    double X,Y;

public:
    // constructeur & destructeur
	Point() { X = Y = -1.; }
	Point( double xi, double yi);
	~Point();

	// fonctions
	void affiche(void)const;
	void attrib_coord(double x, double y);
    void getCoord(double &x,double &y);

};

// Fonctions externes mais qui sont en rapport avec la classe Point
Point* LecPoints(const char* name, int &N);
void afficheCoordPoints(Point* ListPoints, int N);
void CreateFilePoint(const char* name, int nrow, Point *A);
void CreateFileMatPoint(const char* name, int nrow, int ncol, Point **A);

#endif // POINT_H
