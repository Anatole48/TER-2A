#include <fstream>
#include <iostream>
#include "matrix.h"

using namespace Eigen;
using namespace std;


int main()
{
    int tmax = 10;
    double dt = 0.01;
    int n = int (tmax/dt);
    int xmin = 0;
    int  xmax = 1;
    double dx = 0.01;
    int  ymin =0;
    int ymax = 1;
    double dy = 0.01;
    double t = 2;
    int Nx = int(ceil((xmax-xmin)/dx)-1);
    dx = (xmax-xmin)/(Nx+1.);
    int Ny = int(ceil((ymax-ymin)/dy)-1);
    dy = (ymax-ymin)/(Ny+1.);
    SparseMatrix<double> S;
    SparseMatrix<double> A(Nx*Ny,Nx*Ny);
    VectorXd T;
    VectorXd F;
    VectorXd P;
    VectorXd B(Nx*Ny);
    SimplicialLDLT <SparseMatrix<double> > solver;
    SparseMatrix<double> I (Nx*Ny,Nx*Ny);

    vector<Triplet<double>> triplets;
    for (int k =0 ; k<Nx*Ny;k++)
    {
        triplets.push_back({k,k,1});
    }
    I.setFromTriplets(triplets.begin(), triplets.end());
    BuildLaplacianMatrix(xmin ,xmax, ymin,ymax ,dx,dy,dt,S);//on contruit la matrice
    CL(xmin ,xmax, ymin,ymax ,dx,dy,dt,T); // On construit le vecteur contentant les CL
    A = I + S;
    solver.compute(A);
//    cout << A << endl;

    // on contruit le vecteur solution à l'instant intitial
    buildinitial(xmin ,xmax, ymin,ymax,dx,dy,P);

    

    for (int i = 0 ;i<n;i++)
         {
             BuildSource(xmin,xmax,ymin,ymax ,dx ,dy , dt*(i+1),dt,F); //on construit le second membre
             B=F+T+P;

             P = solver.solve(B);

    }
    cout << P << endl;
    ofstream mon_flux;
    string name_file("result.txt");
    mon_flux.open(name_file, ios::out);

    if(mon_flux)
      {
        for (int k = 0 ; k < Nx*Ny  ; k++ )
      {
          int i = k%Nx + 1 ;
          int j = k/Nx + 1;
        mon_flux << i*dx << " " << j*dy<< " " << P[k] << " " << endl;
      }
      }
    else
      {
        cout << "ERREUR: Impossible d’ouvrir le fichier." << endl;
      }
    mon_flux.close();
  

  return 0;

                               
                               
}
