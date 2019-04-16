#include "boundaryCondition.h"

void zeroBound(vectorField &vfield)
{
	for (size_t j = 0; j < rank_y; j++)
	{
		vfield.xval[0][j] = 0;
		vfield.xval[rank_x - 1][j] = 0;

		vfield.yval[0][j] = 0;
		vfield.yval[rank_x - 1][j] = 0;

		vfield.zval[0][j] = 0;
		vfield.zval[rank_x - 1][j] = 0;
	}
}

void zeroBound(scalarField &sfield)
{
	for (size_t j = 0; j < rank_y; j++)
	{
		sfield.val[0][j] = 0;
		sfield.val[rank_x - 1][j] = 0;
	}
}

void periodicBound(vectorField &vfield)
{
	for (size_t i = 0; i < rank_x; i++)
	{
		vfield.xval[i][0] = vfield.xval[i][rank_y - 2];
		vfield.xval[i][rank_y - 1] = vfield.xval[i][1];
		vfield.yval[i][0] = vfield.yval[i][rank_y - 2];
		vfield.yval[i][rank_y - 1] = vfield.yval[i][1];
		vfield.zval[i][0] = vfield.zval[i][rank_y - 2];
		vfield.zval[i][rank_y - 1] = vfield.zval[i][1];   
	}
}

void periodicBound(scalarField &sfield)
{
	for (size_t i = 0; i < rank_x; i++)
	{
		sfield.val[i][0] = sfield.val[i][rank_y - 2];
		sfield.val[i][rank_y - 1] = sfield.val[i][1];
	}
}