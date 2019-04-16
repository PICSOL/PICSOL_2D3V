#pragma once
#include "stdafx.h"
#include "Particle.h"
#include "Vector3D.h"

/*
borisPusher - standard Boris particle pusher in three dimensional velocity space, two dimensional real space
              update V^(n-1/2) to V^(n+1/2), r^(n) to r^(n+1),in which

			  alpha = (q * dt)/(2 * m),   beta = (2 * alpha)/(1 + alpha^2 * B^2)  
*/

//void borisPusher(Particle<double> &part, double alpha, double beta, double delta_t, Vector3D<double> E, Vector3D<double> B);
void borisPusher(Particle<double> &part, double &alpha, double &beta, double &delta_t, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz, int loop);

void borisPusherRewind(Particle<double> &part, double &alpha, double &beta, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz);

/*
borisPusher - modified Boris particle pusher in three dimensional velocity space, two dimensional real space
              update V^(n-1/2) to V^(n+1/2), r^(n) to r^(n+1), in which

			  alpha = (q * dt)/(2 * m),   beta = 1 / (1 + alpha^2 * B^2)  
*/
//void modifiedBorisPusher(Particle<double> &part, double alpha, double beta, double delta_t, Vector3D<double> E, Vector3D<double> B);
void modifiedBorisPusher(Particle<double> &part, double &alpha, double &beta, double &delta_t, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz, int loop);

void modifiedBorisPusherRewind(Particle<double> &part, double &alpha, double &beta, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz);