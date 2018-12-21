#include "Triangle.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip> // pour setw

using namespace std;

int** lectTriangles(const char* name, int &NbTri){
    /* Recuperation du nombre de triangle (NbTri) dans le fichier name et creation de la 
    matrice NT contenant les sommet de chaque triangle (NbTri*3) obtenu dans le fichier name */

    ifstream fichier(name);
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
    /* Pour chaque triangle, recherche les n° de triangles voisins 
           la valeur -1 est affectée s'il n'y a pas de voisin 
    Renvoie la matrive NTV (int NbTri*3) en sortie */

    int** NTV=CreateMat<int>(NbTri,3);

    for(int k=0; k<NbTri; k++){
        NTV[k][0]=trianglevoisin(NbTri,NT[k][1],NT[k][2],k,NT);  // triangle voisin des sommets A2A3
        NTV[k][1]=trianglevoisin(NbTri,NT[k][0],NT[k][2],k,NT);  // triangle voisin des sommets A1A3
        NTV[k][2]=trianglevoisin(NbTri,NT[k][0],NT[k][1],k,NT);  // triangle voisin des sommets A1A2
    }
    return NTV;
}

void CalculCoeff(double &a1, double &b1, double &c1, Point A, Point B){
    // Calcule les coefficents a, b, c de la droite (D): ax+by+c=0 passant par les point A et B

    double Xa, Ya, Xb, Yb;
    A.getCart(Xa,Ya);
    B.getCart(Xb,Yb);
    b1 = Xb - Xa;
    a1 = Ya - Yb;
    c1 = -Ya*b1 - Xa*a1;
}

Point IntersectionDroites(double a1, double a2, double b1, double b2, double c1, double c2){
    // Calcule les coordonnées X et Y d'intersection des deux droites de coeff a,b,c

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
	// Calcule la distance euclidienne entre les Points A et B

	double Xa, Ya, Xb, Yb;
	A.getCart(Xa,Ya);
	B.getCart(Xb,Yb);

	return ( sqrt( (Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb) ) );
}

Point* initOmega(Point* ListPoints, int NbTri, int **NT){
    /* Recherche les Points omega : centre de cercle inscrit de chaque triangle
    Renvoie le vecteur Omega (Point NbTri)*/

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
    /* Recherche les Points Mi : intersection entre le coté d'un triangle et le segment entre 
    les centres de cercles inscrit de triangles voisins.
    Renvoie la matrice NM (Point NbTri*3)*/

    Point** NM=CreateMat<Point>(nbtri,3);

    for (int k=0; k<nbtri; k++){
        double a1, b1, c1;
        double a2, b2, c2;

        for (int i=0; i<3; i++){
            if (NTV[k][i]!=-1){ // Si il y a un triangle voisin
                CalculCoeff(a1,b1,c1,ListPoints[NT[k][(i+1)%3]-1],ListPoints[NT[k][(i+2)%3]-1]); 
                // Coefficients de la droite (D1): a1*x+b1*y+c1=0 qui passe par les sommets A2 et A3
                CalculCoeff(a2,b2,c2,omega[k],omega[NTV[k][i]]);
                // Coefficients de la droite (D2): a2*x+b2*y+c2=0 qui passe par les sommets le centre du cercle 
                // inscrit du triangle k et par le centre du cercle inscrit du triangle voisin de k par les sommets A2 et A3
                NM[k][i]=IntersectionDroites(a1,a2,b1,b2,c1,c2);

            }else { // S'il n'y a pas de triangle voisin
                double X, Y, X1, Y1, X2, Y2;
                ListPoints[NT[k][(i+1)%3]-1].getCart(X1,Y1);
                ListPoints[NT[k][(i+2)%3]-1].getCart(X2,Y2);
                X=fabs(X1-X2)/2+min(X1,X2);
                Y=fabs(Y1-Y2)/2+min(Y1,Y2);
                NM[k][i].attrib_coord(X,Y);
            }
        }
    }

    // calcule les coordonnées barycentriques des points Mi
    CoordBaryMi(NT, ListPoints, NM, nbtri);

    return NM;
}

void CoordBaryMi(int **NT, Point *ListPoints, Point **NM, int NbTri){
    /* Calcule les coordonnées barycentriques des points Mi de chaque triangle */

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

double** ComputeAllCoeff(int NbTri, int NumFonc, int **NT, Point *ListPoints, Point **NM, Point *Omega){
    // NumFonc : numéro de la fonction a interpoler :
    //      - 1 = exponentielle
    //      - 2 = polynomiale

    double** AllCoeffs = new double*[NbTri];

    for (int k = 0; k < NbTri; ++k)
    {
        AllCoeffs[k] = new double[19];  // Vecteur contenant tous les coefficients pour un triangle A1A2A3: ai,bi,ci,di,ei,mi,w

        double p, q, r, a, b, alpha, w1, w2, w3;
        for (int i=0; i<3; i++){
            p=multPointsCart(gradf(ListPoints[NT[k][i]-1],NumFonc),calcVect(ListPoints[NT[k][i]-1],NM[k][(i+1)%3]));  //Initialisation du point pi
            q=multPointsCart(gradf(ListPoints[NT[k][i]-1],NumFonc),calcVect(ListPoints[NT[k][i]-1],NM[k][(i+2)%3]));  //Initialisation du point qi
            r=multPointsCart(gradf(ListPoints[NT[k][i]-1],NumFonc),calcVect(ListPoints[NT[k][i]-1],Omega[k]));  //Initialisation du point ri

            AllCoeffs[k][i]=f(ListPoints[NT[k][i]-1],NumFonc);   // Initialisation des coefficients ai
            AllCoeffs[k][i+3]=AllCoeffs[k][i]+q/2;   //Initialisation des coefficients bi
            AllCoeffs[k][i+6]=AllCoeffs[k][i]+p/2;   //Initialisation des coefficients ci
            AllCoeffs[k][i+9]=AllCoeffs[k][i]+r/2;   //Initialisation des coefficients di
        }

        for (int i=0; i<3; i++){
            switch (i){
            case 0:
                NM[k][0].getBary(b,a,alpha);
                break;
            case 2:
                NM[k][1].getBary(alpha,b,a);
                break;
            case 1:
                NM[k][2].getBary(a,alpha,b);
                break;
            }

            AllCoeffs[k][i+12]=alpha*AllCoeffs[k][9+(i+2)%3]+(1-alpha)*AllCoeffs[k][9+(i+1)%3];
            /*if (k==8){
            cout<<"Ei: "<<AllCoeffs[k][i+12]<<endl;
            cout<<"alpha: "<<alpha<<endl;
            cout<<"1er: "<<alpha*AllCoeffs[k][9+(i+2)%3]<<", 2nd: "<<(1-alpha)*AllCoeffs[k][9+(i+1)%3]<<endl;
            //cout<"dk: "<<AllCoeffs[k][9+(i+2)%3]<< ", dj : "<<AllCoeffs[k][9+(i+1)%3]<<endl;
            }*/
            AllCoeffs[k][i+15]=alpha*AllCoeffs[k][6+(i+2)%3]+(1-alpha)*AllCoeffs[k][3+(i+1)%3];
        }
        Omega[k].getBary(w1,w2,w3);
        AllCoeffs[k][18]=w1*AllCoeffs[k][9]+w2*AllCoeffs[k][10]+w3*AllCoeffs[k][11];
    }

    return AllCoeffs;
}

double foncInterpolant(Point M, Point** MicroTriangles, int l, double *Coeffs){
    // Fonction qui calcule la valeur de L'interpolant à partir du numéro de micro triangle
    // l : numéro du micro-triangle
    // Coeffs : liste des 19 coefficients du triangle T

    int i,j,k;
    bool cas(l%2==0); // True si l = 2i

    if (cas) { i=(l)/2; }
    else { i = (l-1)/2; }
    j = (i+1)%3;
    k = (i+2)%3;

    double Interpolant;
    double w1,w2,w3;

    CartToBary(M,MicroTriangles[l][0],MicroTriangles[l][1],MicroTriangles[l][2]);
    M.getBary(w1,w2,w3);

    if (!cas)
    {
        double mk = Coeffs[15+k];
        double ek = Coeffs[12+k];
        double w = Coeffs[18];
        double bi = Coeffs[3+i];
        double di = Coeffs[9+i];
        double ai = Coeffs[0+i];
        Interpolant = mk*w3*w3 + ek*2*w2*w3 + w*w2*w2 + bi*2*w1*w3 + di*2*w1*w2 + ai*w1*w1;
    }
    else
    {
        double mj = Coeffs[15+j];
        double ej = Coeffs[12+j];
        double w = Coeffs[18];
        double ci = Coeffs[6+i];
        double di = Coeffs[9+i];
        double ai = Coeffs[0+i];
        Interpolant = mj*w3*w3 + ej*2*w2*w3 + w*w2*w2 + ci*2*w1*w3 + di*2*w1*w2 + ai*w1*w1;
    }

    return Interpolant;
}

double evalInterpolant(Point I, Point *ListPoints, int **NT, Point *Omega, Point **NM, double **AllCoeff, Point*** SMT, int NbTri){
    // fonction qui évalue la valeur de l'interpolant au point quelconque I
    // SMT : pointeur des sommets des microtriangle de chaque triangle

    int numTriangle = LocatePointTriangle(I, ListPoints, NT, NbTri);
    int numMicroTriangle = LocatePointMicroTriangle(I, ListPoints, numTriangle, SMT, NT, Omega, NM);

    double val = foncInterpolant(I, SMT[numTriangle], numMicroTriangle, AllCoeff[numTriangle]);

    return val;
}

Point*** ComputeAllSMT(int NbTri, int **NT, Point *ListPoints, Point *Omega, Point **NM){
    // evalue et retourne un pointeur sur les sommets de chaque microtriangle de chaque triangle

    Point*** SMT = new Point**[NbTri];

    for (int k = 0; k < NbTri; ++k)
    {
        SMT[k] = CreateMat<Point>(6,3);
        int j(0);

        for (int i=0; i<3; i++){
            SMT[k][i+j][0]=ListPoints[NT[k][i]-1]; SMT[k][i+j][1]=Omega[k]; SMT[k][i+j][2]=NM[k][(i+1)%3];
            SMT[k][i+1+j][0]=ListPoints[NT[k][i]-1]; SMT[k][i+1+j][1]=Omega[k]; SMT[k][i+1+j][2]=NM[k][(i+2)%3];
            j++;
        }
    }

    return SMT;
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

int LocatePointMicroTriangle(Point A, Point *ListPoints, int k, Point*** SMT, int **NT, Point *Omega, Point **NM){
    //La fonction renvoie le numero du micro triangle où est localise le Point A dans le triangle k

    bool test(false);
    int t(0);
    while(test==false && t<6){
        test=dansMicroTriangle(A,t,SMT[k]);
        t+=1;
    }
    return t-1;
}

void FichierSurface(const char* n_grille , const char* n_fonc, int NumFonc, Point *ListPoint, int NbPts, int **NT, Point *Omega, Point **NM, 
                    int NbTri, double **AllCoeff, Point ***SMT){

    ofstream interp(n_grille);
    ofstream fonc(n_fonc);

    // Calcul des a, b, c et d qui definissent la taille du domaine
    double a, b, c, d;
    double x0, y0, xi, yi;
    ListPoint[0].getCart(x0,y0);
    a=x0; b=x0;
    c=y0; d=y0;
    for(int i=1; i<NbPts; i++){
        ListPoint[i].getCart(xi,yi);
        a=min(a,xi);
        b=max(b,xi);
        c=min(c,yi);
        d=max(d,yi);
    }

    double hx((b-a)/100); //pas entre chaque point de la grille par rapport a l'axe x
    double hy((d-c)/100); //pas entre chaque point de la grille par rapport a l'axe y
    Point grille;
    double x(0.), y(0.);
    for(int i=0; i<100; i++){    //boucle sur l'axe x
        x = (double) (i*hx);
        for(int j=0; j<100; j++){   // boucle sur l'axe y
            y = (double)(j*hy);
            grille.attrib_coord(x,y);

            fonc << f(grille, NumFonc) << " ";
            interp << evalInterpolant(grille,ListPoint,NT,Omega,NM,AllCoeff,SMT,NbTri) << " ";
        }
        fonc << endl;
        interp << endl;
    }


    interp.close();
    fonc.close();
}

///////////////// OUTPUT

void results_ListTriangles(ostream &fichier, int NbTri, int **NT, Point *Omega){
    double X,Y;
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "    Liste des triangles et des points Omega associés" << endl; 
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << setw(5) << "k" << setw(10) << "NT(1,k)" << setw(10) << "NT(2,k)" << setw(10) << "NT(3,k)" << setw(10) << "X_omega" 
            << setw(10) << "Y_omega" << endl;
    for (int k=0; k<NbTri; k++){
        Omega[k].getCart(X,Y);
        fichier << setw(5) << k << setw(10) << NT[k][0] << setw(10) << NT[k][1] << setw(10) << NT[k][2] << setw(10) << X 
                << setw(10) << Y << endl;
    }
}
void results_ListPoints(ostream &fichier, int NbTri, Point **NM){
    double X, Y, X1, X2, X3, Y1, Y2, Y3;
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "    Liste des points Mi pour chaque triangle" << endl; 
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << setw(5) << "k" << setw(10) << "x(M1)" << setw(10) << "y(M1)" << setw(10) << "x(M2)" << setw(10) << "y(M2)" 
            << setw(10) << "x(M3)" << setw(10) << "y(M3)" << endl;
    for (int k=0; k<NbTri; k++){
        NM[k][0].getCart(X1,Y1);
        NM[k][1].getCart(X2,Y2);
        NM[k][2].getCart(X3,Y3);
        fichier << setw(5) << k << setw(10) << X1 << setw(10) << Y1 << setw(10) << X2 << setw(10) << Y2 << setw(10) << X3 
                << setw(10) << Y3 << endl;
    }
}
void results_ValFonc(ostream &fichier, int NumFonc, int NbPts, Point *ListPoints){
    double X,Y;
    fichier << endl << "----------------------- FONCTION N°" << NumFonc << " -------------------------------" << endl << endl;
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "    Liste des N points et valeurs de la fonction" << endl; 
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << setw(5) << "i" << setw(10) << "xi" << setw(10) << "yi" << setw(10) << "f(xi,yi)" << setw(10) << "df/dx" 
            << setw(10) << "df/dy" << endl;
    for (int i=0; i<NbPts; i++){
        ListPoints[i].getCart(X,Y);
        fichier << setw(5) << i << setw(10) << X << setw(10) << Y << setw(10) << f(ListPoints[i],NumFonc) << setw(10) 
                << fpx(ListPoints[i], NumFonc) << setw(10) << fpy(ListPoints[i],NumFonc) << endl;
    }
}
void results_Interpol(ostream &fichier, int NbTri, int NumFonc, Point *ListPoints, int **NT, Point *Omega, Point **NM, 
                      double **AllCoeff, Point*** SMT){
    Point UN(2.5,0.8), DEUX(0.2,1.1), TROIS(2.9,2.5);
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "    Valeur de l'interpolant aux points :" << endl; 
    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "(2.5,0.8) : " << evalInterpolant(UN, ListPoints, NT, Omega, NM, AllCoeff, SMT, NbTri) << endl;
    fichier << "(0.2,1.1) : " << evalInterpolant(DEUX, ListPoints, NT, Omega, NM, AllCoeff, SMT, NbTri) << endl;
    fichier << "(2.9,2.5) : " << evalInterpolant(TROIS, ListPoints, NT, Omega, NM, AllCoeff, SMT, NbTri) << endl;
}
void results_Erreur(ostream &fichier, int NumFonc, Point *ListPoint, int NbPts, int **NT, Point *Omega, Point **NM, 
                    int NbTri, double **AllCoeff, Point ***SMT){

    fichier << "--------------------------------------------------------------------" << endl;
    fichier << "    Valeur minimum et maximum de f-S sur les points de visualisation" << endl; 
    fichier << "--------------------------------------------------------------------" << endl;

    // Calcul des a, b, c et d qui definissent la taille du domaine
    double a, b, c, d;
    double x0, y0, xi, yi;
    ListPoint[0].getCart(x0,y0);
    a=x0; b=x0;
    c=y0; d=y0;
    for(int i=1; i<NbPts; i++){
        ListPoint[i].getCart(xi,yi);
        a=min(a,xi);
        b=max(b,xi);
        c=min(c,yi);
        d=max(d,yi);
    }

    double hx((b-a)/100); //pas entre chaque point de la grille par rapport a l'axe x
    double hy((d-c)/100); //pas entre chaque point de la grille par rapport a l'axe y

    Point grille;
    double x, y, F, S;
    double erreur_max(-100), erreur_min(100);

    for(int i=0; i<100; i++){    //boucle sur l'axe x
        x = (double) (i*hx);
        for(int j=0; j<100; j++){   // boucle sur l'axe y
            y = (double)(j*hy);
            grille.attrib_coord(x,y);

            F = f(grille, NumFonc);
            S = evalInterpolant(grille,ListPoint,NT,Omega,NM,AllCoeff,SMT,NbTri);

            erreur_min = min(erreur_min,(F-S));
            erreur_max = max(erreur_max,(F-S));
        }
    }

    fichier << "erreur min = " << erreur_min << endl;
    fichier << "erreur max = " << erreur_max << endl;
}


