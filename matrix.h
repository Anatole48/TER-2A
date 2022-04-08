#include "EigenLibrary/Eigen/Dense"
#include "EigenLibrary/Eigen/Sparse"
#include "Function.h"


//FOnction contruisant la matrice
void BuildLaplacianMatrix(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy , Eigen ::  SparseMatrix<double>  &S);

//Fonction construisant le terme source
void BuildSource(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy , double t ,Eigen :: VectorXd   &F);
// Fonction contruisant le vecteur solution au temps 0
void buildinitial(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy ,Eigen :: VectorXd   &F);
// FOnction des CL aux bords
void CL(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy ,Eigen :: VectorXd   &F);
