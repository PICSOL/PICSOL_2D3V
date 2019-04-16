#pragma once
#include "stdafx.h"
#include "scalarField.h"
#include "vectorField.h"
#include "Species.h"
#include "Global.h"

/*
linearGather - Specialized function gathers particle information in the cell using linear 
               interploration, also generate particle density field and particle current field
			   and vector field for Et computation.

			   this function is responsible for gathering all needed information about the particles.
*/
void linearGather(Species& species, scalarField& Density, vectorField& Current, vectorField& Convect);

/*
boundEff - modified boundary conditions considering the effective area at the boundary
*/
void boundEff(scalarField &sfield);
void boundEff(vectorField &vfield);