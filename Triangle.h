#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Point.h"
#include "Fonctions_F_G.h"
#include "Matrice.h"
#include <iostream>
#include <fstream>
#include <math.h>

int** lectTriangles(const char* name, int &NbTri);

int trianglevoisin(int nbtri, int S1, int S2, int k, int **NT);

int** initNTV(int NbTri, int **NT);

void CalculCoeff(double &a1, double &b1, double &c1, Point A, Point B);

Point IntersectionDroites(double a1, double a2, double b1, double b2, double c1, double c2);

double DistEucl(Point A, Point B);

Point* initOmega(Point* ListPoints, int NbTri, int **NT);

Point** initNM(int **NT, int **NTV, Point *ListPoints, int nbtri, Point *omega);

void CoordBaryMi(int **NT, Point* ListPoints, Point **NM, int NbTri);

void CreatFileResults(const char* name,int **NT, Point *Omega, Point **NM, Point *ListPoints,int NbTri, int NbPts);

void CartToBary( Point &A, Point S1, Point S2, Point S3);

void BaryToCart( Point &A, Point S1, Point S2, Point S3);

bool dansTriangle(Point& A, int k, int **NT, Point *ListPoints);

int LocatePointTriangle(Point A, double a, double b, double c, double d, Point *ListPoints, int **NT);

double* CoefInterpolation(int k, int **NT, Point *ListPoints, Point **NM, Point *Omega);

#endif // TRIANGLE_H
