#pragma once
/*
Particle - data structure to store position and velocity for one single particle in 2d-3V space
*/

template <class DataType>
class Particle
{
public:
	
	/* variable list */
	DataType x;              // x - position
	DataType y;              // y - position
	DataType vx;             // vx - velocity
	DataType vy;             // vy - velocity
	DataType vz;             // vz - velocity

	/* function list */
	Particle() :x(0),y(0), vx(0),vy(0),vz(0) {}                                                                               // constructor function 
	Particle(DataType x1, DataType y1, DataType vx1, DataType vy1, DataType vz1) :x(x1), y(y1),vx(vx1),vy(vy1),vz(vz1) {}     // alternative constructor function
};

