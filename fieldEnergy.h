#pragma once
#include "vectorField.h"
#include "stdafx.h"
#include "Global.h"

/*
emEnergy - calculate the total electromagnetic field energy, normalized to
           
*/

double emEnergy(vectorField &E_l, vectorField &E_t, vectorField &B);