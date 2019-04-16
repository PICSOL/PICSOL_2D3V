#pragma once
#include "stdafx.h"
#include "Global.h"

/*
tridiagSolver - solving tri-diagonal matrix using Liewellyn Thomas's algorithm, in this case a,b,c are constant arrays
               
			    A(a=1, b, c=1)*x = y, in which A is a constant tri-diagonal matrix, b is the diagnoal element  
*/
void tridiagSolver(double b, vector<double> &x, vector<double> &y);

/*
poissonFFT - 2D poisson solver using fourier transformation method, x is the leading index, dirichlet boundary condition in x direction,
            periodic boundary condition in y direction, using spectrum method to solve poisson equation
			
		    Delta P(x,y) = Q(x,y) in which P(x=0,y) = psi(y), P(x=Lx,y) = phi(y), P(x,0) = P(x, Ly)
*/
void poissonFFT(vector<vector<double>> &P, vector<vector<double>> &Q, vector<double> &LBm_real,
	vector<double> &LBm_imag, vector<double> &RBm_real, vector<double> &RBm_imag);

/*
poissonSOR - 2D poisson solver using SOR iteration method, x is the leading index, dirichlet boundary condition in x direction,
             periodic boundary condition in y direction

             Delta P(x,y) = Q(x,y) in which P(x=0,y) = psi(y), P(x=Lx,y) = phi(y), P(x,0) = P(x, Ly)
*/
int poissonSOR(vector<vector<double>> &P, vector<vector<double>> &Q);


// The FFT method can only be applied to invariable coefficient elliptic equation, however in Darwin-Vlasov model, a variable 
// coefficient is involved, an iterative relaxation method is implemented for solving Helmholtz equation. 

/*
helmholtzSOR - 2D helmholtz solver using SOR iteration method, x is the leading index, dirichlet boundary condition in x direction,
              periodic boundary condition in y direction

              Delta P(x,y) + a(x,y)*P(x,y) = Q(x,y) in which P(x=0,y) = psi(y), P(x=Lx,y) = phi(y), P(x,0) = P(x, Ly)
*/
void helmholtzSOR(vector<vector<double>> &P, vector<vector<double>> &Q, vector<vector<double>> &a);

