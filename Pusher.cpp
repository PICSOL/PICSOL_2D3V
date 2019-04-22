#include "Pusher.h"

void borisPusher(Particle<double> &part, double &alpha, double &beta, double &delta_t, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz, int loop)
{
	double electContrib_x = alpha * Ex;
	double electContrib_y = alpha * Ey;
	double electContrib_z = alpha * Ez;

	double vx_minus, vy_minus, vz_minus, vx_prime, vy_prime, vz_prime;

	for (int i = 0; i < loop; i++)
	{
		vx_minus = part.vx + electContrib_x;
		vy_minus = part.vy + electContrib_y;
		vz_minus = part.vz + electContrib_z;

		vx_prime = vx_minus + alpha * (vy_minus*Bz - vz_minus * By);
		vy_prime = vy_minus + alpha * (vz_minus*Bx - vx_minus * Bz);
		vz_prime = vz_minus + alpha * (vx_minus*By - vy_minus * Bx);

		part.vx = vx_minus + electContrib_x + beta * (vy_prime*Bz - vz_prime * By);
		part.vy = vy_minus + electContrib_y + beta * (vz_prime*Bx - vx_prime * Bz);
		part.vz = vz_minus + electContrib_z + beta * (vx_prime*By - vy_prime * Bx);

		part.x += part.vx * delta_t;
		part.y += part.vy * delta_t;
	}
}

void borisPusherRewind(Particle<double> &part, double &alpha, double &beta, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz)
{
	double electContrib_x = alpha * Ex;
	double electContrib_y = alpha * Ey;
	double electContrib_z = alpha * Ez;

	double vx_minus = part.vx + electContrib_x;
	double vy_minus = part.vy + electContrib_y;
	double vz_minus = part.vz + electContrib_z;

	double vx_prime = vx_minus + alpha * (vy_minus*Bz - vz_minus * By);
	double vy_prime = vy_minus + alpha * (vz_minus*Bx - vx_minus * Bz);
	double vz_prime = vz_minus + alpha * (vx_minus*By - vy_minus * Bx);

	part.vx = vx_minus + electContrib_x + beta * (vy_prime*Bz - vz_prime * By);
	part.vy = vy_minus + electContrib_y + beta * (vz_prime*Bx - vx_prime * Bz);
	part.vz = vz_minus + electContrib_z + beta * (vx_prime*By - vy_prime * Bx);
}

void modifiedBorisPusher(Particle<double> &part, double &alpha, double &beta, double &delta_t, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz, int loop)
{
	double electContrib_x = alpha * Ex;
	double electContrib_y = alpha * Ey;
	double electContrib_z = alpha * Ez;

	double vx_minus, vy_minus, vz_minus, vx_prime, vy_prime, vz_prime, temp;

	for (int i = 0; i < loop; i++)
	{
		vx_minus = part.vx + alpha * electContrib_x;
		vy_minus = part.vy + alpha * electContrib_y;
		vz_minus = part.vz + alpha * electContrib_z;

		temp = vx_minus * Bx + vy_minus * By + vz_minus * Bz;

		vx_prime = beta * (vx_minus + alpha * ((vy_minus*Bz - vz_minus * By) + alpha * temp * Bx));
		vy_prime = beta * (vy_minus + alpha * ((vz_minus*Bx - vx_minus * Bz) + alpha * temp * By));
		vz_prime = beta * (vz_minus + alpha * ((vx_minus*By - vy_minus * Bx) + alpha * temp * Bz));

		part.vx = 2 * vx_prime - part.vx;
		part.vy = 2 * vy_prime - part.vy;
		part.vz = 2 * vz_prime - part.vz;

		part.x += part.vx * delta_t;
		part.y += part.vy * delta_t;
	}
}

void modifiedBorisPusherRewind(Particle<double> &part, double &alpha, double &beta, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz)
{
	double vx_minus = part.vx + alpha * Ex;
	double vy_minus = part.vy + alpha * Ey;
	double vz_minus = part.vz + alpha * Ez;

	double temp = vx_minus * Bx + vy_minus * By + vz_minus * Bz;

	double vx_prime = beta * (vx_minus + alpha * ((vy_minus*Bz - vz_minus * By) + alpha * temp * Bx));
	double vy_prime = beta * (vy_minus + alpha * ((vz_minus*Bx - vx_minus * Bz) + alpha * temp * By));
	double vz_prime = beta * (vz_minus + alpha * ((vx_minus*By - vy_minus * Bx) + alpha * temp * Bz));

	part.vx = 2 * vx_prime - part.vx;
	part.vy = 2 * vy_prime - part.vy;
	part.vz = 2 * vz_prime - part.vz;
}

void pusher(Species& species, vectorField& E, vectorField& B, double dt, int loop)
{	
	/* boris pusher coefficient */
	double alpha = 0.5 * species.specific_charge * dt * lambda;  // lambda is the normalization factor
	double alpha2 = pow(alpha, 2);
	double ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4;
	double Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, beta;
	int i, j, j_plus;

	/* record the removing particle index */
	vector<int> removeIndex;

	int particle_size = species.part.size();
	/* OpenMP multithreading area for particle iteration */
    #pragma omp parallel for private(ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4,Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, beta,i, j, j_plus)
	for (int n = 0; n < particle_size; n++)
	{
	    ratio_x = (species.part[n].x - lbound) * idx;
		i = (int)floor(ratio_x);
		prop_x1 = ratio_x - i;

		ratio_y = (species.part[n].y - dbound) * idy;
		j = (int)floor(ratio_y);
		j_plus = (j == (rank_y - 1)) ? 0 : (j + 1);
	    prop_y1 = ratio_y - j;

		/* 2D linear interploration coefficients */
		a4 = prop_x1 * prop_y1;
		a2 = prop_y1 - a4;
		a3 = prop_x1 - a4;
		a1 = 1.0 - a4 - a3 - a2;
		
		/* linear interploration */
		Exeff = E.xval[i][j] * a1 + E.xval[i][j_plus] * a2 + E.xval[i + 1][j] * a3 + E.xval[i + 1][j_plus] * a4;
		Eyeff = E.yval[i][j] * a1 + E.yval[i][j_plus] * a2 + E.yval[i + 1][j] * a3 + E.yval[i + 1][j_plus] * a4;
		Ezeff = E.zval[i][j] * a1 + E.zval[i][j_plus] * a2 + E.zval[i + 1][j] * a3 + E.zval[i + 1][j_plus] * a4;
		Bxeff = B.xval[i][j] * a1 + B.xval[i][j_plus] * a2 + B.xval[i + 1][j] * a3 + B.xval[i + 1][j_plus] * a4;
		Byeff = B.yval[i][j] * a1 + B.yval[i][j_plus] * a2 + B.yval[i + 1][j] * a3 + B.yval[i + 1][j_plus] * a4;
		Bzeff = B.zval[i][j] * a1 + B.zval[i][j_plus] * a2 + B.zval[i + 1][j] * a3 + B.zval[i + 1][j_plus] * a4;
		
		/* advance particle using Boris method */
		beta = (2 * alpha) / (1 + alpha2 * (Bxeff * Bxeff + Byeff * Byeff + Bzeff * Bzeff));
		borisPusher(species.part[n], alpha, beta, dt, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, loop);

		/* absorbing wall in x direction, periodic in y direction */
		if (species.part[n].x < lbound || species.part[n].x > rbound)
		{
			removeIndex.push_back(n);
		}
		else if (species.part[n].y < dbound)
		{
			species.part[n].y = Ly - fmod(dbound - species.part[n].y, Ly) + dbound;
		}
		else if (species.part[n].y > ubound)
		{
			species.part[n].y = fmod(species.part[n].y - ubound, Ly) + dbound;
		}
	}

	int removeSize = removeIndex.size();
	/* remove the particles that strike the conducting wall */
	for (int i = removeSize - 1; i >= 0; i--)
	{
		/* delete particle from the back */
		swap(species.part[removeIndex[i]],species.part.back());
		species.part.pop_back();
	}
}

void leapFrogRewind(Species& species, vectorField& E, vectorField& B, double rewindTime)
{
	/* boris pusher coefficient */
	double alpha = 0.5 * species.specific_charge * rewindTime * lambda;  // lambda is the normalization factor
	double alpha2 = pow(alpha, 2);
	double ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4;
	double Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, beta;
	int i, j, j_plus;

	/* record the removing particle index */
	vector<int> removeIndex;

	int particle_size = species.part.size();
	/* OpenMP multithreading area for particle iteration */
    #pragma omp parallel for private(ratio_x, ratio_y, prop_x1, prop_y1, a1, a2, a3, a4,Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, beta,i, j, j_plus)
	for (int n = 0; n < particle_size; n++)
	{
		ratio_x = (species.part[n].x - lbound) * idx;
		i = (int)floor(ratio_x);
		prop_x1 = ratio_x - i;

		ratio_y = (species.part[n].y - dbound) * idy;
		j = (int)floor(ratio_y);
		j_plus = (j == (rank_y - 1)) ? 0 : (j + 1);
		prop_y1 = ratio_y - j;

		/* 2D linear interploration coefficients */
		a4 = prop_x1 * prop_y1;
		a2 = prop_y1 - a4;
		a3 = prop_x1 - a4;
		a1 = 1.0 - a4 - a3 - a2;

		/* linear interploration */
		Exeff = E.xval[i][j] * a1 + E.xval[i][j_plus] * a2 + E.xval[i + 1][j] * a3 + E.xval[i + 1][j_plus] * a4;
		Eyeff = E.yval[i][j] * a1 + E.yval[i][j_plus] * a2 + E.yval[i + 1][j] * a3 + E.yval[i + 1][j_plus] * a4;
		Ezeff = E.zval[i][j] * a1 + E.zval[i][j_plus] * a2 + E.zval[i + 1][j] * a3 + E.zval[i + 1][j_plus] * a4;
		Bxeff = B.xval[i][j] * a1 + B.xval[i][j_plus] * a2 + B.xval[i + 1][j] * a3 + B.xval[i + 1][j_plus] * a4;
		Byeff = B.yval[i][j] * a1 + B.yval[i][j_plus] * a2 + B.yval[i + 1][j] * a3 + B.yval[i + 1][j_plus] * a4;
		Bzeff = B.zval[i][j] * a1 + B.zval[i][j_plus] * a2 + B.zval[i + 1][j] * a3 + B.zval[i + 1][j_plus] * a4;

		/* advance particle using Boris method */
		beta = (2 * alpha) / (1 + alpha2 * (Bxeff * Bxeff + Byeff * Byeff + Bzeff * Bzeff));
		borisPusherRewind(species.part[n], alpha, beta, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff);

		/* absorbing wall in x direction, periodic in y direction */
		if (species.part[n].x < lbound || species.part[n].x > rbound)
		{
			removeIndex.push_back(n);
		}
		else if (species.part[n].y < dbound)
		{
			species.part[n].y = Ly - fmod(dbound - species.part[n].y, Ly) + dbound;
		}
		else if (species.part[n].y > ubound)
		{
			species.part[n].y = fmod(species.part[n].y - ubound, Ly) + dbound;
		}
	}

	int removeSize = removeIndex.size();
	/* remove the particles that strike the conducting wall */
	for (int i = removeSize - 1; i >= 0; i--)
	{
		/* delete particle from the back */
		swap(species.part[removeIndex[i]], species.part.back());
		species.part.pop_back();
	}
}