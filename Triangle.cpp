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
    A.getCart(Xa,Ya);
    B.getCart(Xb,Yb);
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
	A.getCart(Xa,Ya);
	B.getCart(Xb,Yb);

	return ( sqrt( (Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb) ) );
}

Point* initOmega(Point* ListPoints, int NbTri, int **NT){
    // Initialisation de la matrice Omega de Nbtri Points

    Point *Omega=new Point[NbTri];
    double w1,w2,w3,sumw;
    double X1,X2,X3,Y1,Y2,Y3;
    double X, Y;

    for (int k=0; k<NbTri; k++){

        // Recuperation des coordonnees des sommets A1, A2 et A3 du triangle k
        Point A1=ListPoints[NT[k][0]-1];
        Point A2=ListPoints[NT[k][1]-1];
        Point A3=ListPoints[NT[k][2]-1];

        A1.getCart(X1,Y1);
        A2.getCart(X2,Y2);
        A3.getCart(X3,Y3);

        // Coordonnees barycentrique de omega du triangle k
        w1 = DistEucl(A2,A3);
        w2 = DistEucl(A1,A3);
        w3 = DistEucl(A1,A2);
        sumw = w1+w2+w3;
        Omega[k].attrib_bary((w1/sumw),(w2/sumw),(w3/sumw));

        // Coordonnnees cartesiennes de omega du triangle k
        X = (w1*X1 + w2*X2 + w3*X3) / sumw;
        Y = (w1*Y1 + w2*Y2 + w3*Y3) / sumw;
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
                ListPoints[NT[k][(i+1)%3]-1].getCart(X1,Y1);
                ListPoints[NT[k][(i+2)%3]-1].getCart(X2,Y2);
                X=fabs(X1-X2)/2+min(X1,X2);
                Y=fabs(Y1-Y2)/2+min(Y1,Y2);
                NM[k][i].attrib_coord(X,Y);
            }
        }
    }
    return NM;
}

void CoordBaryMi(int **NT, Point *ListPoints, Point **NM, int NbTri){
    // Mise en memoire des coordonnees barycentriques des points M1, M2 et M3

    Point A1,A2,A3,M1,M2,M3;
    double alpha1, alpha2, alpha3;
    for (int i = 0; i < NbTri; ++i){
        // Recuperation des coordonnees des sommets A1, A2 et A3 du triangle i
        A1 = ListPoints[NT[i][0]-1];
        A2 = ListPoints[NT[i][1]-1];
        A3 = ListPoints[NT[i][2]-1];

         // Recuperation des coordonnées des 3 Mi du triangle i
        M1 = NM[i][0];
        M2 = NM[i][1];
        M3 = NM[i][2];

        alpha1 = DistEucl(A2,M1) / DistEucl(A2,A3);
        alpha2 = DistEucl(A1,M3) / DistEucl(A1,A2);
        alpha3 = DistEucl(A3,M2) / DistEucl(A3,A1);

         // Attribution des coordonnées barycentriques aux points Mi
        NM[i][0].attrib_bary(0, 1-alpha1, alpha1);
        NM[i][1].attrib_bary(alpha2, 0, 1-alpha2);
        NM[i][2].attrib_bary(1-alpha3, alpha3, 0);
     }
}

void CreatFileResults(const char* name,int **NT, Point *Omega, Point **NM, Point *ListPoints,int NbTri, int NbPts){
    // Creation du fichier PS.RES
    double X, Y, X1, X2, X3, Y1, Y2, Y3;
    ofstream fichier(name);
    for (int k=0; k<NbTri; k++){

        Omega[k].getCart(X,Y);
        fichier<<k<<" "<<NT[k][0]<<" "<<NT[k][1]<<" "<<NT[k][2]<<" "<<X<<" "<<Y<<endl;
    }
    for (int k=0; k<NbTri; k++){
        NM[k][0].getCart(X1,Y1);
        NM[k][1].getCart(X2,Y2);
        NM[k][2].getCart(X3,Y3);
        fichier<<k<<" "<<X1<<" "<<Y1<<" "<<X2<<" "<<Y2<<" "<<X3<<" "<<Y3<<endl;
    }
    for (int i=0; i<NbPts; i++){
        ListPoints[i].getCart(X,Y);
        fichier<<i<<" "<<X<<" "<<Y<<" "<<f(ListPoints[i])<<" "<<fpx(ListPoints[i])<<" "<<fpy(ListPoints[i])<<endl;
    }
    fichier.close();
}

void CartToBary(Point& A, Point S1, Point S2, Point S3){
    // calcule les coordonnées barycentriques du point A par rapport au triangle de sommets S1,S2,S3

    double w1,w2,w3;
    double x1,x2,x3,y1,y2,y3,X,Y;

    // recupere les coordonnées cartésiennes
    S1.getCart(x1, y1);
    S2.getCart(x2, y2);
    S3.getCart(x3, y3);
    A.getCart(X, Y);

    w1 = ( (y2-y3)*(X-x3)+(x3-x2)*(Y-y3) ) / ( (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3) );
    w2 = ( (y3-y1)*(X-x3)+(x1-x3)*(Y-y3) ) / ( (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3) );
    w3 = 1 - w1 - w2;

    // attribution des coordonnées barycentriques au point A
    A.attrib_bary(w1, w2, w3);
}

void BaryToCart( Point &A, Point S1, Point S2, Point S3){
    // calcule les coordonnées cartesiennes d'un point A par rapport au triangle de sommet S1,S2,S3
     double w1,w2,w3;
    double x1,x2,x3,y1,y2,y3;

    S1.getCart(x1, y1);
    S2.getCart(x2, y2);
    S3.getCart(x3, y3);
    A.getBary(w1,w2,w3);
    double X = w1*x1 + w2*x2 + w3*x3;
    double Y = w1*y1 + w2*y2 + w3*x3;
     A.attrib_coord(X,Y);
}

bool dansTriangle(Point& A, int k, int **NT, Point *ListPoints){
    // Renvoie un booleen si le Point A est dans le triangle k

    CartToBary(A,ListPoints[NT[k][0]-1],ListPoints[NT[k][1]-1],ListPoints[NT[k][2]-1]);
    double w1,w2,w3;
    A.getBary(w1,w2,w3);
    if(w1<0 || w2<0 || w3<0){
        return false;
    }else{
        return true;
    }
}

int LocatePointTriangle(Point A, Point *ListPoints, int **NT, int nbtri){
    //La fonction renvoie le numero du triangle où est localise le Point A dans le domaine D=[a,b]*[c,d]

    int k(0);   //numero du triangle contenant le Point A. On demmarre du triangle 0 pour faire la recherche
    bool test(false);
    while(test==false && k<nbtri){
        test=dansTriangle(A,k,NT,ListPoints);
        k+=1;
    }
    return k-1;
}

double* CoefInterpolation(int k, int **NT, Point *ListPoints, Point **NM, Point *Omega){
    /* Renvoie un vecteur de double qui contient tous les coefficients pour 1 triangle donne en argument
    k= la reference du triangle que l'on souhaite etudier
    NT= matrice qui contient les sommets relatifs a chaque sommet
    ListPoints= vecteur de Point qui comporte les coordonnees cartesiennes de chaque triangle
    NM= matrice de Point qui comporte les coordonnees cartesiennes de M1, M2 et M3 pour chaque triangle
    Omega= vecteur de Point qui comporte les coordonnees cartesiennes de omega pour chaque triangle
    */

    double *coefInter=new double[19];  // Vecteur contenant tous les coefficients pour un triangle A1A2A3: ai,bi,ci,di,ei,mi,w

    double p, q, r, a, b, alpha, w1, w2, w3;
    for (int i=0; i<3; i++){
        p=multPointsCart(gradf(ListPoints[NT[k][i]]),calcVect(ListPoints[NT[k][i]],NM[k][(i+1)%3]));  //Initialisation du point pi
        q=multPointsCart(gradf(ListPoints[NT[k][i]]),calcVect(ListPoints[NT[k][i]],NM[k][(i+2)%3]));  //Initialisation du point qi
        r=multPointsCart(gradf(ListPoints[NT[k][i]]),calcVect(ListPoints[NT[k][i]],Omega[k]));  //Initialisation du point ri

        coefInter[i]=f(ListPoints[NT[k][i]]);   // Initialisation des coefficients ai
        coefInter[i+3]=coefInter[i]+q/2;   //Initialisation des coefficients bi
        coefInter[i+6]=coefInter[i]+p/2;   //Initialisation des coefficients ci
        coefInter[i+9]=coefInter[i]+r/2;   //Initialisation des coefficients di

        switch (i){
        case 0:
            NM[k][0].getBary(b,a,alpha);
        case 1:
            NM[k][1].getBary(alpha,b,a);
        case 2:
            NM[k][2].getBary(a,alpha,b);
        }
        Omega[k].getBary(w1,w2,w3);

        coefInter[i+12]=alpha*coefInter[9+(i+2)%3]+(1-alpha)*coefInter[9+(i+1)%3];
        coefInter[i+15]=alpha*coefInter[6+(i+2)%3]+(1-alpha)*coefInter[3+(i+1)%3];
    }
    coefInter[18]=w1*coefInter[9]+w2*coefInter[10]+w3*coefInter[11];

    return coefInter;
}

Point** MicroTriangle(int k, int **NT, Point *ListPoints, Point *Omega, Point **NM){
    // Renvoie une matrice de Point (6*3) qui renvoie les sommets de chaque microtriangle dans un triangle

    Point **NMT=CreateMat<Point>(6,3);
    int j(0);

    // A VERIFIER CETTE PARTIE CAR C'EST ICI QU'IL Y A UN BEUG
    for (int i=0; i<3; i++){
        NMT[i+j][0]=ListPoints[NT[k][i]-1]; NMT[i+j][1]=Omega[k]; NMT[i+j][2]=NM[k][(i+1)%3];
        NMT[i+1+j][0]=ListPoints[NT[k][i]-1]; NMT[i+1+j][1]=NM[k][(i+2)%3]; NMT[i+1+j][2]=Omega[k];
        j++;
    }

    return NMT;
}

bool dansMicroTriangle(Point& A, int t, Point **NMT){
    // Renvoie un booleen si le Point A est dans le triangle k

    CartToBary(A,NMT[t][0],NMT[t][1],NMT[t][2]);
    double w1,w2,w3;
    A.getBary(w1,w2,w3);
    if(w1<0 || w2<0 || w3<0){
        return false;
    }else{
        return true;
    }
}

int LocatePointMicroTriangle(Point A, Point *ListPoints, int k, int **NT, Point *Omega, Point **NM){
    //La fonction renvoie le numero du micro triangle où est localise le Point A dans le triangle k

    Point **NMT=MicroTriangle(k,NT,ListPoints,Omega,NM);
    bool test(false);
    int t(0);
    while(test==false && t<6){
        test=dansMicroTriangle(A,t,NMT);
        t+=1;
    }
    return t-1;
}


