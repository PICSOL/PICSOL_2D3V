#include "fieldEnergy.h"

double emEnergy(vectorField &E_total, vectorField &B)
{
	double total_energy = 0;

	for (size_t i = 1; i < rank_x - 1; i++)
	{
		for (size_t j = 0; j < rank_y; j++)
		{
			total_energy += (E_total.xval[i][j] * E_total.xval[i][j]
				+ E_total.yval[i][j] * E_total.yval[i][j]
				+ E_total.zval[i][j] * E_total.zval[i][j])*lambda
				+ B.xval[i][j] * B.xval[i][j]
				+ B.yval[i][j] * B.yval[i][j]
				+ B.zval[i][j] * B.zval[i][j];
		}   
	}

	/* dirichlet boundary only have half effective area */
	for (size_t j = 0; j < rank_y; j++)
	{
		total_energy += 0.5 * (E_total.xval[0][j] * E_total.xval[0][j]
			+ E_total.yval[0][j] * E_total.yval[0][j]
			+ E_total.zval[0][j] * E_total.zval[0][j])*lambda
			+ B.xval[0][j] * B.xval[0][j]
			+ B.yval[0][j] * B.yval[0][j]
			+ B.zval[0][j] * B.zval[0][j];
	}

	for (size_t j = 0; j < rank_y; j++)
	{
		total_energy += 0.5 * (E_total.xval[rank_x - 1][j] * E_total.xval[rank_x - 1][j]
			+ E_total.yval[rank_x - 1][j] * E_total.yval[rank_x - 1][j]
			+ E_total.zval[rank_x - 1][j] * E_total.zval[rank_x - 1][j])*lambda
			+ B.xval[rank_x - 1][j] * B.xval[rank_x - 1][j]
			+ B.yval[rank_x - 1][j] * B.yval[rank_x - 1][j]
			+ B.zval[rank_x - 1][j] * B.zval[rank_x - 1][j];
	}

	return total_energy;
}