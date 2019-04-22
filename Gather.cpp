#include "Gather.h"

void linearGather(Species& species, scalarField& Density, vectorField& Current, vectorField& Convect)
{   
	/* field value reset to zero */
	Density.clear();
	Current.clear();
	Convect.clear();

	double ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4;
	double vx, vy, vz, t1, t2, b1, b2, b3, b4;
	int i, j, j_plus;
	
	/* constant parameter */
	static const double weight = (double)rank_x * rank_y / (double)initParticle;
	int particle_size = species.part.size();

	/* OpenMP multithreading area for particle iteration */
    #pragma omp parallel for private(ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4, vx, vy, vz, t1, t2, b1, b2, b3, b4, i, j, j_plus)
	for (int n = 0; n < particle_size; n++)
	{
		ratio_x = (species.part[n].x - lbound) * idx;
		i = (int)floor(ratio_x);
		prop_x1 = ratio_x - i;
	

		/* periodic in y direction */
		ratio_y = (species.part[n].y - dbound) * idy;
		j = (int)floor(ratio_y);
		j_plus = (j == (rank_y - 1)) ? 0 : (j + 1);
		prop_y1 = ratio_y - j;
		

		/* 2D linear interploration coefficients */
		a4 = prop_x1 * prop_y1;
		a2 = prop_y1 - a4;
		a3 = prop_x1 - a4;
		a1 = 1 - a4 - a3 - a2;

		/* density computation */
        #pragma omp atomic
		Density.val[i][j] += a1;
        #pragma omp atomic
		Density.val[i][j_plus] += a2;
        #pragma omp atomic
		Density.val[i + 1][j] += a3;
        #pragma omp atomic
		Density.val[i + 1][j_plus] += a4;

		/* current computation */
		vx = species.part[n].vx;
		vy = species.part[n].vy;
		vz = species.part[n].vz;

        #pragma omp atomic
		Current.xval[i][j] += vx * a1;
        #pragma omp atomic
		Current.yval[i][j] += vy * a1;
        #pragma omp atomic
		Current.zval[i][j] += vz * a1;

        #pragma omp atomic
		Current.xval[i][j_plus] += vx * a2;
        #pragma omp atomic
		Current.yval[i][j_plus] += vy * a2;
        #pragma omp atomic
		Current.zval[i][j_plus] += vz * a2;
        #pragma omp atomic
		Current.xval[i + 1][j] += vx * a3;
        #pragma omp atomic
		Current.yval[i + 1][j] += vy * a3;
        #pragma omp atomic	
		Current.zval[i + 1][j] += vz * a3;
        #pragma omp atomic
		Current.xval[i + 1][j_plus] += vx * a4;
        #pragma omp atomic
		Current.yval[i + 1][j_plus] += vy * a4;
        #pragma omp atomic
		Current.zval[i + 1][j_plus] += vz * a4;

		/* convect computation */
		t1 = idx * vx;
		t2 = idy * vy;

		b1 = -t1 - t2;
		b2 = -t1 + t2;
		b3 = -b2;
		b4 = -b1;

        #pragma omp atomic
		Convect.xval[i][j] += vx * b1;
        #pragma omp atomic
		Convect.yval[i][j] += vy * b1;
        #pragma omp atomic
		Convect.zval[i][j] += vz * b1;
        #pragma omp atomic
		Convect.xval[i][j_plus] += vx * b2;
        #pragma omp atomic
		Convect.yval[i][j_plus] += vy * b2;
        #pragma omp atomic
		Convect.zval[i][j_plus] += vz * b2;
        #pragma omp atomic
		Convect.xval[i + 1][j] += vx * b3;
        #pragma omp atomic
		Convect.yval[i + 1][j] += vy * b3;
        #pragma omp atomic
		Convect.zval[i + 1][j] += vz * b3;
        #pragma omp atomic
		Convect.xval[i + 1][j_plus] += vx * b4;
        #pragma omp atomic
		Convect.yval[i + 1][j_plus] += vy * b4;
        #pragma omp atomic
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