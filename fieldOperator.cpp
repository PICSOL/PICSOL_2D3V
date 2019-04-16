#include "fieldOperator.h"

void gradient(vectorField &grad, scalarField &sfield, double coef)
{
	/* coefficient table */
	double coefx = 0.5 * idx * coef;
	double coefy = 0.5 * idy * coef;

	for (size_t i = 1; i < rank_x - 1; i++)
	{
		for (size_t j = 1; j < rank_y - 1; j++)
		{
			grad.xval[i][j] = coefx * (sfield.val[i + 1][j] - sfield.val[i - 1][j]);
			grad.yval[i][j] = coefy * (sfield.val[i][j + 1] - sfield.val[i][j - 1]);
		}

		grad.xval[i][0] = coefx * (sfield.val[i + 1][0] - sfield.val[i - 1][0]);
		grad.xval[i][rank_y - 1] = coefx * (sfield.val[i + 1][rank_y - 1] - sfield.val[i - 1][rank_y - 1]);
		grad.yval[i][0] = coefy * (sfield.val[i][1] - sfield.val[i][rank_y - 1]);
		grad.yval[i][rank_y - 1] = coefy * (sfield.val[i][0] - sfield.val[i][rank_y - 2]);
	}
}

void divergence(scalarField &div, vectorField &vfield, double coef)
{
	/* coefficient table */
	double coefx = 0.5 * idx * coef;
	double coefy = 0.5 * idy * coef;

	for (size_t i = 1; i < rank_x - 1; i++)
	{
		for (size_t j = 1; j < rank_y - 1; j++)
		{
			div.val[i][j] = coefx * (vfield.xval[i + 1][j] - vfield.xval[i - 1][j])
				          + coefy * (vfield.yval[i][j + 1] - vfield.xval[i][j - 1]);
		}
		
		div.val[i][0] = coefx * (vfield.xval[i + 1][0] - vfield.xval[i - 1][0])
			          + coefy * (vfield.yval[i][1] - vfield.xval[i][rank_y - 1]);

		div.val[i][rank_y -1] = coefx * (vfield.xval[i + 1][rank_y - 1] - vfield.xval[i - 1][rank_y - 1])
			                  + coefy * (vfield.yval[i][0] - vfield.xval[i][rank_y - 2]);
	}
}

void rotation(vectorField &rot, vectorField &vfield, double coef) 
{
	/* coefficient table */
	double coefx1 = -0.5 * idx * coef;
	double coefx2 = 0.5 * idx * coef;
	double coefy1 = -0.5 * idy * coef;
	double coefy2 = 0.5 * idy * coef;

	for (size_t i = 1; i < rank_x - 1; i++)
	{
		for (size_t j = 1; j < rank_y - 1; j++)
		{
			rot.xval[i][j] = coefy2 * (vfield.zval[i][j + 1] - vfield.zval[i][j - 1]);
			rot.yval[i][j] = coefx1 * (vfield.zval[i + 1][j] - vfield.zval[i - 1][j]);
			rot.zval[i][j] = coefx2 * (vfield.yval[i + 1][j] - vfield.yval[i - 1][j])
				           + coefy1 * (vfield.xval[i][j + 1] - vfield.xval[i][j - 1]);
		}

		rot.xval[i][0] = coefy2 * (vfield.zval[i][1] - vfield.zval[i][rank_y - 1]);
		rot.xval[i][rank_y - 1] = coefy2 * (vfield.zval[i][0] - vfield.zval[i][rank_y - 2]);

		rot.yval[i][0] = coefx1 * (vfield.zval[i + 1][0] - vfield.zval[i - 1][0]);
		rot.yval[i][rank_y - 1] = coefx1 * (vfield.zval[i + 1][rank_y - 1] - vfield.zval[i - 1][rank_y - 1]);
	
		rot.zval[i][0] = coefx2 * (vfield.yval[i + 1][0] - vfield.yval[i - 1][0])
			           + coefy1 * (vfield.xval[i][1] - vfield.xval[i][rank_y - 1]);

		rot.zval[i][rank_y - 1] = coefx2 * (vfield.yval[i + 1][rank_y - 1] - vfield.yval[i - 1][rank_y - 1])
			                    + coefy1 * (vfield.xval[i][0] - vfield.xval[i][rank_y - 2]);
	}
}