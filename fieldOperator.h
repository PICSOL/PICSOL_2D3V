#pragma once
#include "Global.h"
#include "scalarField.h"
#include "vectorField.h"

/*
*  Implement for three basic field operators: grad, div, rot, 
*  those operators are valid only under the circumstances that z direction is 
*  neglectable, which implies 'partial_z = 0' and boundary condition for x direction
*  must be dirichlet boundary condition, boundary condition for y direction must be 
*  periodic boundary condition. coef is the constant coefficient.
*/

/*
grad - gradient operator for scalarfield
*/
void gradient(vectorField &grad, scalarField &sfield, double coef);

/*
div - divergence operator for vectorfield
*/
void divergence(scalarField &div, vectorField &vfield, double coef);

/*
rot - rotation operator for vectorfield
*/
void rotation(vectorField &rot, vectorField &vfield, double coef);

