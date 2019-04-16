#pragma once
#include "Global.h"
#include "scalarField.h"
#include "vectorField.h"


/*
smooth2D - smooth function for scalar field and vector field using binomial smoothing method
           specialized for dirichlet B.C.in x direction, periodic B.C.in y direction
*/
void smooth2D(scalarField &field);
void smooth2D(vectorField &field);


/*
smoothParticleField - smooth function for particle related field using binomial smoothing method
                      specialized for dirichlet B.C.in x direction, periodic B.C.in y direction
*/
void smoothParticleField(scalarField &electDen,
                         scalarField &ionDen,
                         vectorField &electConvect,
                         vectorField &ionConvect,
                         vectorField &electCurrent,
                         vectorField &ionCurrent);