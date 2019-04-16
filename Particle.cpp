#include "Particle.h"

/* sort funtion for Particle based on the first variable */
template <class DataType>
static bool particleSort(const Particle<DataType>& a1, const Particle<DataType>& a2)
{
	return a1.x <= a2.x;
}