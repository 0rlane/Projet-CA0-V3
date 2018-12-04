#ifndef FONCTION_F_G_H
#define FONCTION_F_G_H

#include "Point.h"
#include "Matrice.h"
#include <iostream>
#include <fstream>
#include <math.h>

double f(Point A);

double fpx(Point A);

double fpy(Point A);

double g(Point A);

double gpx(Point A);

double gpy(Point A);

double** initF(int NbPts, Point* ListPoints);

double** initG(int NbPts, Point* ListPoints);


#endif // FONCTION_F_G_H
