#pragma once
#include "stdafx.h"
#include "vectorField.h"
#include "Species.h"
#include "particlePusher.h"
#include "Global.h"

/*
pusher - function to update the position and velocity of one species 
         based on the field quantities on the grid-point and boundary condition,
		 this function must be coherent with "gather" function
*/
void pusher(Species& species, vectorField& E, vectorField& B, double dt, int loop);

/*
leapFrogRewind - function to rewind particle velocity to form a staggered leap frog format
*/
void leapFrogRewind(Species& species, vectorField& E, vectorField& B, double rewindTime);
