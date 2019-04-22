#include "stdafx.h"
#include "scalarField.h"
#include "vectorField.h"
#include "Particle.h"
#include "Species.h"
#include "Param.h"
#include "fieldUpdater.h"
#include "fieldEnergy.h"
#include "Pusher.h"
#include "Smooth.h"
#include "Gather.h"
#include "auxiliary.h"
#include "boundaryCondition.h"
#include "Global.h"

/********************************************************************************************************* 
 PICSOL is a 2d-3v fully kinetic parallel electromagnetic PIC code developed by Huaxiang Zhang since 2018 Mar.28th,
 which is primarily designed to simulate low-recycling rate SOL(Scrape-off Layer) plasma biasing in tokamaks. 
 PICSOL is built based on 2D Vlasov-Darwin model, which is capable of discribing plasma's kinetic behavior without 
 electromagnetic wave. This code is open source to anyone who intends to use it. 
 The latest version can be downloaded from https://github.com/ZhangPlasma/PICSOL/
 Version Number: 1.03
**********************************************************************************************************/

int main() {   
/*********************************************************************************************************/
/*********************************************************************************************************/
/* Input laoding */

    /*initialize global variables */
	initialGlobal(param);

	/* create data storage file folder based on system time */
	string directory = recordFile(param.PATH);
	
	/* create log file for simulation run */
	ofstream log_file = logFile(directory);
	param.saveInput(directory);
	
	/* RNG seeding */
	srand((unsigned int)time(NULL));

/***********************************************************************************************************/
/***********************************************************************************************************/
/* Initialization */
	
	/* particle initialization */
	Species electron(param.ELECTMASS, param.ELECTCHARGE, param.ETEMP, param.ETEMP, param.ETEMP, param.DELTA_Te);
	Species ion(param.IONMASS, param.IONCHARGE, param.ITEMP, param.ITEMP, param.ITEMP, param.DELTA_Ti);
	
	/* maxwell distribution initialization */	
	electron.initialMaxwell(param);
	ion.initialMaxwell(param);
    
	/* initial plasma kinetic energy */
	double totalElectKinetic = electron.kinetic();
	double totalIonKinetic = ion.kinetic();
	
	/* log_file record initial state */
	recordInitial(electron, ion, log_file);

	/* electromagnetic field declaration */
	vectorField E_l(param);
	vectorField E_t(param);
	vectorField E_total(param);
	vectorField B(param);
	scalarField Phi(param);

	/* particle related field declaration */
	scalarField electDen(param);
	scalarField ionDen(param);
	vectorField electConvect(param);
	vectorField ionConvect(param);
	vectorField electCurrent(param);
	vectorField ionCurrent(param);
    vectorField J(param);
	scalarField rho(param);
	vectorField convect(param);

	/* particle field initialization */
	linearGather(electron,electDen,electCurrent,electConvect);
	linearGather(ion, ionDen, ionCurrent, ionConvect);
	smoothParticleField(electDen, ionDen, electConvect, ionConvect, electCurrent, ionCurrent);
	rhoComp(electDen, ionDen, rho);
	JComp(electCurrent, ionCurrent, J);
	convectComp(electConvect, ionConvect, convect);
	
	/* electromagnetic field initialization */
	phiComp(rho, Phi);
	ElComp(Phi, E_l);
	//BComp(J, B); 
	B.add(B_0);
	//EtComp(rho, J, convect, E_l, B, E_t);
	EtotalComp(E_l, E_t, E_total);

	/* use electron advancing time as the minimal time step */
	int nu = (int)floor(param.DELTA_T / electron.dt);
	double field_dt = nu * electron.dt;

	/* leapfrog rewinding */
	leapFrogRewind(electron, E_total, B, -0.5 * electron.dt);
	leapFrogRewind(ion, E_total, B, -0.5 * field_dt);

	/* save initial state */
	recordLog(param.INITIAL_PARTICLE, 0, electron, ion, totalElectKinetic, totalIonKinetic, log_file);
	recordData(0, directory, electron, ion, electCurrent, electDen, ionCurrent, ionDen, Phi, rho, E_l, E_t, B);

/*******************************************************************************************************/
/*******************************************************************************************************/
/* MAIN LOOP */

	clock_t start = clock();

	for (int main_iter = 1; main_iter <= param.MAX_LOOP; main_iter++)
	{	
		/* sub-stepped particle pusher */
		pusher(electron, E_total, B, electron.dt, nu);
		pusher(ion, E_total, B, field_dt, 1);
         
		/* particle field update */
		linearGather(electron, electDen, electCurrent, electConvect);
		linearGather(ion, ionDen, ionCurrent, ionConvect);
		smoothParticleField(electDen, ionDen, electConvect, ionConvect, electCurrent, ionCurrent);
		rhoComp(electDen, ionDen, rho);
		JComp(electCurrent, ionCurrent, J);
		convectComp(electConvect, ionConvect, convect);

		/* electromagnetic field update */
		phiComp(rho, Phi);
		ElComp(Phi, E_l);
		//BComp(J, B);
		//EtComp(rho, J, convect, E_l, B, E_t);
		EtotalComp(E_l, E_t, E_total);

		/* save data */
		if (main_iter % param.RECORD_CHECK == 0)
		{   
			recordLog(param.INITIAL_PARTICLE, main_iter, electron, ion, totalElectKinetic, totalIonKinetic, log_file);
			recordData(main_iter, directory, electron, ion, electCurrent, electDen, ionCurrent, ionDen, Phi, rho, E_l, E_t, B);
		}
	}
	clock_t end = clock();
	log_file << "Main Loop is complete, for " << param.MAX_LOOP << " steps, time consuming: " << 1.0 * (end - start) / (double)(CLOCKS_PER_SEC) << "s" << endl;
	return 0;
}

