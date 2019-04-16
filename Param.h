#pragma once
#include "stdafx.h"
/*
Param - data structure to store simulation parameter
*/

class Param
{
public:

	/* physical constant */
	const double  EPS_0 = 8.85418782e-12;    // F/m, vacuum permittivity
	const double  MU_0 = 1.256637061e-6;     // N/A^2, vacuum permeability 
	const double  K = 1.38065e-23;			 // J/K, Boltzmann constant
	const double  ME = 9.10938215e-31;   	 // kg, electron mass
	const double  QE = 1.602176565e-19;		 // C, elementary charge
	const double  AMU = 1.660538921e-27;     // kg, atomic mass unit
	const double  PI = 3.14159265358979;     // circumference ratio
	const double  e = 2.718281828459045;     // base of natural logarithm

	/* data storage path */
	string PATH;

	/* simulation flag */
	int  SHAPE_FACTOR;        // 0 for NGP, 1 for CIC, 2 for second order spline

	/* simulation constant */	
	int  MAX_ITERATION;       // maximum iteration for Gauss-Seidel method
	double  MIN_ERROR;        // minimum error tolerance for iteration method    

	int  GRID_NUM_X;          // grid point number in x direction
	int  GRID_NUM_Y;          // grid point number in y direction
	int  MAX_LOOP;            // maximum cycle number
	int  RECORD_CHECK;        // iteration number to record data
	int  INITIAL_PARTICLE;    // initial particle number
	int  INITIAL_PARTICLE_X;  // initial particle number in x direction
	int  INITIAL_PARTICLE_Y;  // initial particle number in y direction
	double  DELTA_T;          // minimum time step for field update, normlized to DELTA_Te
	double  DELTA_Te;         // minimum time step for electron update, normalized to inverse of plasma frequency 1/omega_{pe}
	double  DELTA_Ti;         // minimum time step for ion update, normalized to DELTA_Te
	double  DELTA_X;          // minimum space step in x direction, normalized to electron Debye length Lambda_{De}
	double  DELTA_Y;          // minimum space step in y direction, normalized to electron Debye length Lambda_{De}
	double  LBOUND;           // simulation domain left bound in x direction
	double  RBOUND;           // simulation domain right bound in x direction
	double  UBOUND;           // simulation domain upper bound in y direction
	double  DBOUND;           // simulation domain lower bound in y direction
	double  XLENGTH;          // simulation domain length in x direction
	double  YLENGTH;          // simulation domain length in y direction 
	double  BX;               // constant background magnetic field in x direction, in Tesla
	double  BY;               // constant background magnetic field in y direction, in Tesla 
	double  REFLUX_SURFACE;   // refluxing surface position

	/* plasma parameters */
	double  DENSITY;            // characteristic plasma density 
	double  ETEMP;              // electron temperature, normalized to electron temperature
	double  ITEMP;              // ion temperature, normalized to electron temperature
	double  ELECTCHARGE;        // electron charge, normalized to elementary charge
	double  IONCHARGE;          // ion charge, normalized to elementary charge
	double  ELECTMASS;          // electron mass, normalized to electron mass
	double  IONMASS;            // ion mass (deuterium), normalized to electron mass

	/* divertor biasing parameters */
	vector<double> LB;          // voltage distribution on the Left bound in x direction
	vector<double> RB;          // voltage distribution on the Right bound in x direction 
	double V_BIAS;              // max biasing voltage
	double Y_0;                 // center of biasing 
	double L_0;                 // biasing length 
	double L_1;                 // length of the slope


    /* function list */
	Param();                                                     // loading simulation constant from given directory
    void copyFile(string& input, string& ouput, string& dir);    // copy file as a binary file
	void saveInput(string& destname);                            // save input data 
};