#include "matrix.h"
#include <iostream>
using namespace std;
using namespace Eigen;
void BuildLaplacianMatrix(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy , Eigen ::  SparseMatrix<double>  &S)
{
    int Nx = int(ceil((xmax-xmin)/dx)-1);
    dx = (xmax-xmin)/(Nx+1.);
    int Ny = int(ceil((ymax-ymin)/dy)-1);
    dy = (ymax-ymin)/(Ny+1.);
    S.resize(Nx*Ny,Nx*Ny);
    vector<Triplet<double>> triplets;
    for (int k=0; k<Nx*Ny; k++)
    {
        int i = k%Nx + 1 ;
        int j = k/Nx + 1;
        double alpha_i_plus_1demi= 1;
        double alpha_j_plus_1demi= 1;
        double alpha_i_moins_1demi= 1;
        double alpha_j_moins_1demi= 1;
        double Sij = 1;

//        cout << i<< j<< endl;

        triplets.push_back({k,k,(alpha_i_plus_1demi+alpha_i_moins_1demi)/((pow(dx,2))*Sij)+ (alpha_j_plus_1demi+alpha_j_moins_1demi)/((pow(dy,2))*Sij)});  // on place la diagonale


        if (i > 1)     // Le sommet n'est pas près d'un bord gauche
        {
            triplets.push_back({k,k-1,-alpha_i_moins_1demi/((pow(dx,2))*Sij)});
            
        }


        if (i < Nx) {  // Le sommet n'est pas près d'un bord droit
            triplets.push_back({k,k+1, -alpha_i_plus_1demi/((pow(dx,2))*Sij)});
        }

        // Le sommet n'est pas près d'un bord haut

        if (j < Ny)
        {
            triplets.push_back({k,k+Nx,-alpha_j_plus_1demi/((pow(dy,2))*Sij)});// on place les b de droite
        }

        // Le sommet n'est pas près d'un bas
        if ( j>1)
        {
            triplets.push_back({k,k-Nx,-alpha_j_moins_1demi/((pow(dy,2))*Sij)});// on place les b de gauche
        }
    }
//    cout << S.size()<< endl;
    S.setFromTriplets(triplets.begin(), triplets.end());
    


}


void BuildSource(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy , double t ,Eigen :: VectorXd   &F)
{
    int Nx = int(ceil((xmax-xmin)/dx)-1);
    dx = (xmax-xmin)/(Nx+1.);
    int Ny = int(ceil((ymax-ymin)/dy)-1);
    dy = (ymax-ymin)/(Ny+1.);
    F.resize(Nx*Ny);
    for (int k = 0 ; k<Nx*Ny;k++)
    {
        int i = k%Nx + 1;
        int j = k/Nx + 1;
        F(k) = f_instat(i*dx,j*dy,t);
    }
    
}
void buildinitial(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy ,Eigen :: VectorXd   &F)
{
    int Nx = int(ceil((xmax-xmin)/dx)-1);
    dx = (xmax-xmin)/(Nx+1.);
    int Ny = int(ceil((ymax-ymin)/dy)-1);
    dy = (ymax-ymin)/(Ny+1.);
    F.resize(Nx*Ny);
    for (int k = 0 ; k<Nx*Ny;k++)
    {
        int i = k%Nx + 1;
        int j = k/Nx + 1;
        F(k) = f(i*dx,j*dy);
    }
    
}

void CL(int xmin, int xmax , int ymin , int ymax ,double &dx , double &dy ,Eigen :: VectorXd   &F)
{
    int Nx = int(ceil((xmax-xmin)/dx)-1);
    dx = (xmax-xmin)/(Nx+1.);
    int Ny = int(ceil((ymax-ymin)/dy)-1);
    dy = (ymax-ymin)/(Ny+1.);
    F.resize(Nx*Ny);
    for (int k = 0 ;k<Nx*Ny;k++)
    {
        double b=0;
        double c=0;
        double d=0;
        double e=0;
        double f = 0;
        int i = k%Nx+1;
        int j = k/Nx +1;
        double alpha_i_plus_1demi= 1;
        double alpha_j_plus_1demi= 1;
        double alpha_i_moins_1demi= 1;
        double alpha_j_moins_1demi= 1;
        double Sij = 1;


        if (i==1) // on est sur un bord gauche
        {
            b = ((alpha_i_moins_1demi)/(dx*dx*Sij)) * h((i-1)*dx,j*dy);
        }
        if (i==Nx) // on est sur un bord droit
        {
            c = ((alpha_i_plus_1demi)/(dx*dx*Sij)) * g((i+1)*dx,j*dy);
        }
        if ( j == 1 ) // on est  sur un bord bas
        {
//            d  = ((alpha_j_moins_1demi)/(dy*dy*Sij)) * h(i*dx,(j-1)*dy);
            d=0;
        }
        if (j==Ny)  // on est sur un bord haut
        {
//            e = ((alpha_j_plus_1demi)/(dy*dy*Sij)) * h(i*dx,(j+1)*dy);
            e=0;
        }
        F[k]= b+c+d+e;
        
        
        
    }
}

