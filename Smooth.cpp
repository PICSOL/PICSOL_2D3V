#include "Smooth.h"

void smooth2D(scalarField &field)
{
    #pragma omp parallel for
	for (int i = 1; i < (int)rank_x - 1; i++)
	{
		/* periodic boundary in y direction */
		field.val[i][0] = 0.125 * (field.val[i][0] * 4 +
			field.val[i][1] +
			field.val[i][rank_y - 1] +
			field.val[i - 1][0] +
			field.val[i + 1][0]);
		
		for (size_t j = 1; j < rank_y - 1; j++)
		{
			field.val[i][j] = 0.125 * (field.val[i][j] * 4 +
				field.val[i][j - 1] +
				field.val[i][j + 1] +
				field.val[i - 1][j] +
				field.val[i + 1][j]);
		}

		field.val[i][rank_y - 1] = 0.125 * (field.val[i][rank_y - 1] * 4 +
			field.val[i][0] +
			field.val[i][rank_y - 2] +
			field.val[i - 1][rank_y - 1] +
			field.val[i + 1][rank_y - 1]);
	}
}

void smooth2D(vectorField &field)
{
    #pragma omp parallel for
	for (int i = 1; i < (int)rank_x - 1; i++)
	{
		/* periodic boundary in y direction */
		field.xval[i][0] = 0.125 * (field.xval[i][0] * 4 +
			field.xval[i][1] +
			field.xval[i][rank_y - 1] +
			field.xval[i - 1][0] +
			field.xval[i + 1][0]);

		field.yval[i][0] = 0.125 * (field.yval[i][0] * 4 +
			field.yval[i][1] +
			field.yval[i][rank_y - 1] +
			field.yval[i - 1][0] +
			field.yval[i + 1][0]);

		field.zval[i][0] = 0.125 * (field.zval[i][0] * 4 +
			field.zval[i][1] +
			field.zval[i][rank_y - 1] +
			field.zval[i - 1][0] +
			field.zval[i + 1][0]);
		
		for (size_t j = 1; j < rank_y - 1; j++)
		{
			field.xval[i][j] = 0.125 * (field.xval[i][j] * 4 +
				field.xval[i][j - 1] +
				field.xval[i][j + 1] +
				field.xval[i - 1][j] +
				field.xval[i + 1][j]);

			field.yval[i][j] = 0.125 * (field.yval[i][j] * 4 +
				field.yval[i][j - 1] +
				field.yval[i][j + 1] +
				field.yval[i - 1][j] +
				field.yval[i + 1][j]);

			field.zval[i][j] = 0.125 * (field.zval[i][j] * 4 +
				field.zval[i][j - 1] +
				field.zval[i][j + 1] +
				field.zval[i - 1][j] +
				field.zval[i + 1][j]);
		}

		field.xval[i][rank_y - 1] = 0.125 * (field.xval[i][rank_y - 1] * 4 +
			field.xval[i][0] +
			field.xval[i][rank_y - 2] +
			field.xval[i - 1][rank_y - 1] +
			field.xval[i + 1][rank_y - 1]);

		field.yval[i][rank_y - 1] = 0.125 * (field.yval[i][rank_y - 1] * 4 +
			field.yval[i][0] +
			field.yval[i][rank_y - 2] +
			field.yval[i - 1][rank_y - 1] +
			field.yval[i + 1][rank_y - 1]);

		field.zval[i][rank_y - 1] = 0.125 * (field.zval[i][rank_y - 1] * 4 +
			field.zval[i][0] +
			field.zval[i][rank_y - 2] +
			field.zval[i - 1][rank_y - 1] +
			field.zval[i + 1][rank_y - 1]);
	}
}

void smoothParticleField(scalarField &electDen, scalarField &ionDen, 
	                     vectorField &electConvect, vectorField &ionConvect,
	                     vectorField &electCurrent, vectorField &ionCurrent)
{
	smooth2D(electDen);
	smooth2D(ionDen);
	smooth2D(electConvect);
	smooth2D(ionConvect);
	smooth2D(electCurrent);
	smooth2D(ionCurrent);
}