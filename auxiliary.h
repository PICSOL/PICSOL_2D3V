#pragma once
#include "stdafx.h"
#include "Species.h"
#include "scalarField.h"
#include "vectorField.h"

/*
recordFile - record file directory
*/
string recordFile(string path);

/*
logFile - create log file 
*/
ofstream logFile(string directory);

/*
recordInitial - record particle information in log file
*/
void recordInitial(Species& electron, Species& ion, ofstream& log_file);

/*
recordLog - record field and particle information in log file
*/
void recordLog(double INITIAL_PARTICLE, int main_iter, Species& electron, Species& ion, double totalElectKinetic, double totalIonKinetic, ofstream& log_file);

/*
recordData - record field and particle data in data file 
*/
void recordData(int main_iter, string directory, Species& electron, Species& ion, vectorField& electCurrent,
	scalarField& electDen, vectorField& ionCurrent, scalarField& ionDen, scalarField& phi, scalarField& rho, vectorField& E_l,
	vectorField& E_t,vectorField& B );
