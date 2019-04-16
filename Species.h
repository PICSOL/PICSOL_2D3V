#pragma once
#include "stdafx.h"
#include "Particle.h"
#include "Param.h"

/*
Species - data structure to store vector and basic info for particular particle species
*/
class Species
{
public:

	/* variable list */
	double mass;                          // particle's mass, normalized to electron mass me
	double charge;                        // particle's charge, normalized to elementary charge e
	double specific_charge;               // specific charge, normalized to e/me
	double dt;                            // minimum time advance step 
	double tx;                            // particle's average temperature in x direction, normalized to electron temperature
	double ty;                            // particle's average temperature in y direction, normalized to electron temperature
	double tz;                            // particle's average temperature in z direction, normalized to electron temperature
	vector<Particle<double>> part;        // particle vector to store particle information

	/* function list */
	Species(double mass, double charge, double tx, double ty, double tz, double dt); // constructor function
 	void initialMaxwell(Param &param);                                               // initialize Maxwell distributed particles
	void addParticle(double x, double y, double vx, double vy, double vz);           // add a particle with given velocity and position
	void remParticle(size_t index);                                                  // remove a particle with given index number 
	double kinetic();        	                                                     // return particle's total kinetic energy
	void saveParticle(string directory);                                             // save particle information in given directory
	void saveGlobal(string directory);                                               // save global particle number and kinetic energy                                                         // return the total kinetic energy for this species
};	

/*
generateMaxwell - generate maxwell distridution random number N(mu,sigma) based on Box-Muller Method
*/
double generateMaxwell(double mu, double sigma);

/*
generateMaxwellFlux - generate maxwell flux distridution random number abs(v) * N(sigma) based on inverse transform of sampling
*/
double generateMaxwellFlux(double sigma);

/*
generateMaxwellBeamFlux - generate maxwell with beam flux distridution random number based on inverse transform of sampling
*/
double generateMaxwellBeamFlux(double alpha, double vt, double vb);