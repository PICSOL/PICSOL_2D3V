#pragma once
#include "stdafx.h"
#include "Param.h"
#include "Vector3D.h"

/*
vectorField - data structure to store vector field value using vector container
*/

class vectorField
{
public:

	/* variable list */
	int x_grid_num;                        // number of grid point in x direction
	int y_grid_num;                        // number of grid point in y direction
	double xbound1;                        // left bound of grid region
	double xbound2;                        // right bound of grid region 
	double ybound1;                        // lower bound of grid region
	double ybound2;                        // upper bound of grid region 
	double x_inc;                          // displacement between each grid point in x direction
	double y_inc;                          // displacement between each grid point in y direction
	vector<vector<double>> xval;           // x component using vector container for Field value storage, x is the leading direction, [i][j]'s upper bound is [x_grid_num - 1][y_grid_num - 1]
	vector<vector<double>> yval;           // y component
	vector<vector<double>> zval;           // z component

    /* funtion list */
	vectorField(Param &param);                                           // field initialization using given parameter
	void equal(Vector3D<double> vec);                                    // constant field 
	void multiply(double constant);                                      // grid value multiply by a constant
	void add(Vector3D<double> vec);                                      // add a constant field    
	void sum(vectorField &a, vectorField &b);                            // sum of two vector field
	void substract(vectorField &a);                                      // substraction of two vector field
	void clear();                                                        // reset grid value to zero
	void saveInfo(string directory1,string directory2,string directory3);// save gird value in the given directory
	virtual ~vectorField();                                              // overwrite distructor function to avoid memory leak
};
