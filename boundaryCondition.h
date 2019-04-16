#pragma once
#include "Global.h"
#include "vectorField.h"
#include "scalarField.h"

/*
zeroBound - adjust boundary value for dirichlet boundary condition in x direction 
*/

void zeroBound(vectorField &vfield);
void zeroBound(scalarField &sfield);

/*
periodicBound - periodic boundary in y direction using redundant gird point method
*/

void periodicBound(vectorField &vfield);
void periodicBound(scalarField &sfield);
