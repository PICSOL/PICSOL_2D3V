#pragma once
#include "Param.h"
#include "FFT.h"
#include "Vector3D.h"

/*
 simulation parameters determined by the user input file
*/

/*
Global variables declaration
*/

/* simulation parameters class */
extern Param param;

/* simulation constant */
extern const size_t rank_x;       // number of grid points in x direction
extern const size_t rank_y;       // number of grid points in y direction
extern const double lbound;       // position of the left bound
extern const double rbound;       // position of the right bound
extern const double ubound;       // position of the upper bound
extern const double dbound;       // position of the lower bound
extern const double dx;           // space step in x direction
extern const double idx;          // 1 / dx
extern const double dx2;          // dx^2
extern const double idx2;         // 1 / dx^2
extern const double dy;           // space step in y direction
extern const double idy;          // 1 / dy
extern const double dy2;          // dy^2
extern const double idy2;         // 1 / dy^2
extern const double dxy2;         // -2.0 * (1.0 / dy2 + 1.0 / dx2)
extern const double idxy2;        // 1 / (-2.0 * (1.0 / dy2 + 1.0 / dx2))
extern const double Lx;           // simulation region's length in x direction
extern const double Ly;           // simulation region's length in y direction
extern const double twoPi;        // 2 * Pi
extern const size_t maxIteration; // maximum iteration number for SOR
extern const double minError;     // error tolerance for SOR
extern const double ioncharge;    // ion charge normlized to electron
extern const double electcharge;  // electron charge
extern const int initParticle;    // number of the initial particle 

/* dimensionless number from equation normalization */
extern const double lambda;    // (V_{te} / c)^2

/* normalization constant for electromagnetic field */
extern const double normEl;     // sqrt(T_e * n_e / eps_0) 
extern const double normEt;     // (V_{te} / c)^2 * sqrt(T_e * n_e / eps_0) 
extern const double normB;      // (V_{te} / c^2) * sqrt(T_e * n_e / eps_0)

/* constant background magnetic field */
extern const Vector3D<double> B_0;

/*
omegaTable for the SOR chebyshev acceleration
*/
extern vector<double> omegaTable;

/*
real-valued FFT of phi,psi on the static boundary condition
*/
extern vector<double> phi_m_real;
extern vector<double> phi_m_imag;
extern vector<double> psi_m_real;
extern vector<double> psi_m_imag;

/*
real-valued FFT of zeros on static doundary condition
*/
extern vector<double> FFTbound;

/*
coefficient table diagElement for tridiagSolver in FFT
*/
extern vector<double> diagElement;

/*
coefficient table in Real-valued FFT
*/
extern vector<double> ArTable;
extern vector<double> AiTable;
extern vector<double> BrTable;
extern vector<double> BiTable;

/*
coefficient table in Imaginary-valued FFT
*/
extern vector<double> sinTable;
extern vector<double> cosTable;

/*
initialization for all global tables and variables
*/
extern void initialGlobal();