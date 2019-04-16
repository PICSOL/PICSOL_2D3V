#include "Gather.h"

void linearGather(Species& species, scalarField& Density, vectorField& Current, vectorField& Convect)
{
	/* constant parameter */
	static const double weight = 1.0 * rank_x * rank_y / (double)initParticle;
	
	/* field value reset to zero */
	Density.clear();
	Current.clear();
	Convect.clear();

	/* particle gathering process, CIC(Cloud In Cell) method implement */
	for (size_t n = 0; n < species.part.size(); n++)
	{
		double ratio_x = (species.part[n].x - lbound) * idx;
        int i = (int)floor(ratio_x);
		double prop_x1 = ratio_x - i;
		double prop_x2 = 1.0 - prop_x1;
		
		/* periodic in y direction */
		double ratio_y = (species.part[n].y - dbound) * idy;
		int j = (int)floor(ratio_y);
		int j_plus = (j == (rank_y - 1)) ? 0 : (j + 1);
		double prop_y1 = ratio_y - j;
		double prop_y2 = 1.0 - prop_y1;

		/* 2D linear interploration coefficients */
		double a1 = prop_x2 * prop_y2;
		double a2 = prop_x2 * prop_y1;
		double a3 = prop_x1 * prop_y2;
		double a4 = prop_x1 * prop_y1;

		/* density computation */
		Density.val[i][j] += a1;
		Density.val[i][j_plus] += a2;
		Density.val[i + 1][j] += a3;
		Density.val[i + 1][j_plus] += a4;
        
		/* current computation */
		double vx = species.part[n].vx;
		double vy = species.part[n].vy;
		double vz = species.part[n].vz;

		Current.xval[i][j] += vx * a1;
		Current.yval[i][j] += vy * a1;
		Current.zval[i][j] += vz * a1;

		Current.xval[i][j_plus] += vx * a2;
		Current.yval[i][j_plus] += vy * a2;
		Current.zval[i][j_plus] += vz * a2;

		Current.xval[i + 1][j] += vx * a3;
		Current.yval[i + 1][j] += vy * a3;
		Current.zval[i + 1][j] += vz * a3;

		Current.xval[i + 1][j_plus] += vx * a4;
		Current.yval[i + 1][j_plus] += vy * a4;
		Current.zval[i + 1][j_plus] += vz * a4;

		/* convect computation */
		double t1 = idx * vx;
		double t2 = idy * vy;

		double b1 = - t1 - t2;
		double b2 = - t1 + t2;
		double b3 = - b2;
		double b4 = - b1;

		Convect.xval[i][j] += vx * b1;
		Convect.yval[i][j] += vy * b1;
		Convect.zval[i][j] += vz * b1;

		Convect.xval[i][j_plus] += vx * b2;
		Convect.yval[i][j_plus] += vy * b2;
		Convect.zval[i][j_plus] += vz * b2;

		Convect.xval[i + 1][j] += vx * b3;
		Convect.yval[i + 1][j] += vy * b3;
		Convect.zval[i + 1][j] += vz * b3;

		Convect.xval[i + 1][j_plus] += vx * b4;
		Convect.yval[i + 1][j_plus] += vy * b4;
		Convect.zval[i + 1][j_plus] += vz * b4;
	}
	
	/* normaliztion for physical variables */
	Density.multiply(weight);
	Current.multiply(species.charge * weight);
	Convect.multiply(species.charge * weight);

	/* modifying boundary values due to the effective area */
	boundEff(Density);
	boundEff(Current);
	boundEff(Convect);
}

void boundEff(scalarField &sfield)
{
	/* grid points at dirichlet boundary only have half effective area */
	for (size_t j = 0; j < rank_y; j++)
	{
		sfield.val[0][j] *= 2.0;
		sfield.val[rank_x - 1][j] *= 2.0;
	}
}

void boundEff(vectorField &vfield)
{
	/* grid points at dirichlet boundary only have half effective area */
	for (size_t j = 0; j < rank_y; j++)
	{
		vfield.xval[0][j] *= 2.0;
		vfield.xval[rank_x - 1][j] *= 2.0;
		vfield.yval[0][j] *= 2.0;
		vfield.yval[rank_x - 1][j] *= 2.0;
		vfield.zval[0][j] *= 2.0;
		vfield.zval[rank_x - 1][j] *= 2.0;
	}
}