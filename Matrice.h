#ifndef MATRICE_H
#define MATRICE_H

#include <fstream>
#include <iostream>

template<typename Type> Type** CreateMat(int nrow, int ncol){
    // Creation d'une matrice de taille nrow*ncol

	Type** Mat = new Type*[nrow];
    for (int i = 0; i < nrow; ++i){
    	Mat[i] = new Type[ncol];
    }
    return Mat;
}

template<typename Type> Type AfficheMat(int nrow, int ncol, Type **A){
  // Affichage d'une matrice nrow*ncol

  for (int i=0; i<nrow; i++){
    for(int j=0; j<ncol;j++){
      std::cout<<A[i][j]<<" ";
    }
    std::cout<<std::endl;;
  }
}

template<typename Type> Type FreeMat(Type **A, int nrow){
    // Desalocaation memoire de la matrice A

    for (int i=0; i<nrow; i++) {
            delete [] (A[i]);
    }
    delete[] A;
}

template<typename Type> Type CreateFileMat(const char* name, int nrow, int ncol, Type **A){
    // Creation d'un fichier qui contient la matrice A

	std::ofstream fichier(name);
    for (int i = 0; i < nrow; ++i)
    {
    	for (int j = 0; j < ncol; ++j)
    	{
    		fichier << A[i][j] << " ";
    	}
        fichier << std::endl;
    }
    fichier.close();
}

template<typename Type> Type multVect(Type *A, Type *B, int sA, int sB){
    // Multipliaction du vecteur A de taille sA avec le vecteur B de taille sB
    double C(0);
    if (sA == sB){
        for(int i=0; i<sA; i++){
            C+=A[i]*B[i];
        }
    }else{
        std::cout<<"Erreur dans la mutiplication de vecteur";
    }
    return C;
}

#endif // MATRICE_H

