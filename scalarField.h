#pragma once
#include "stdafx.h"
#include "Param.h"
/*
scalarField - data structure to store scalar field value using vector container
*/

class scalarField
{
public:

	/* variable list */
    int x_grid_num;                // number of grid point in x direction
	int y_grid_num;                // number of grid point in y direction
	double xbound1;                // left bound of grid region
	double xbound2;                // right bound of grid region 
	double ybound1;                // lower bound of grid region
	double ybound2;                // upper bound of grid region 
	double x_inc;                  // displacement between each grid point in x direction
	double y_inc;                  // displacement between each grid point in y direction
	vector<vector<double>> val;    // using vector container for Field value storage, x is the leading direction, [i][j]'s upper bound is [x_grid_num - 1][y_grid_num - 1]

	/* funtion list */
	scalarField(Param &param);                                           // field initialization using given parameter
	void equal(double constant);                                         // constant field 
	void multiply(double constant);                                      // grid value multiply by a constant
	void sum(scalarField& a, scalarField& b);                            // sum of two scalar field
	void clear();                                                        // reset grid value to zero
	void saveInfo(string directory);                                     // save gird value in the given directory
	virtual ~scalarField();                                              // overwrite distructor function to avoid memory leak
};

