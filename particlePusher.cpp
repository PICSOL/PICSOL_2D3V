#include "particlePusher.h"

void borisPusher(Particle<double> &part, double &alpha, double &beta, double &delta_t, double &Ex, double &Ey, double &Ez,
	double &Bx, double &By, double &Bz, int loop)
{
	double electContrib_x = alpha * Ex;
	double electContrib_y = alpha * Ey;
	double electContrib_z = alpha * Ez;
	double vx_minus, vy_minus, vz_minus;
	double vx_prime, vy_prime, vz_prime;

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
	double vx_minus, vy_minus, vz_minus;
	double vx_prime, vy_prime, vz_prime, temp;

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