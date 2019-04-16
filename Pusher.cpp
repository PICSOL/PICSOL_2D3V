#include "Pusher.h"

void pusher(Species& species, vectorField& E, vectorField& B, double dt, int loop)
{
	/* effective electromagnetic field */
	double Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff;
	
	/* boris pusher coefficient */
	double alpha = 0.5 * species.specific_charge * dt * lambda;  // lambda is the normalization factor
	double alpha2 = pow(alpha, 2);

	for (size_t n = 0; n < species.part.size(); n++)
	{
		double ratio_x = (species.part[n].x - lbound) * idx;
		int i = (int)floor(ratio_x);
		double prop_x1 = ratio_x - i;
		double prop_x2 = 1.0 - prop_x1;

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
		
		/* linear interploration */
		Exeff = E.xval[i][j] * a1 + E.xval[i][j_plus] * a2 + E.xval[i + 1][j] * a3 + E.xval[i + 1][j_plus] * a4;
		Eyeff = E.yval[i][j] * a1 + E.yval[i][j_plus] * a2 + E.yval[i + 1][j] * a3 + E.yval[i + 1][j_plus] * a4;
		Ezeff = E.zval[i][j] * a1 + E.zval[i][j_plus] * a2 + E.zval[i + 1][j] * a3 + E.zval[i + 1][j_plus] * a4;
		Bxeff = B.xval[i][j] * a1 + B.xval[i][j_plus] * a2 + B.xval[i + 1][j] * a3 + B.xval[i + 1][j_plus] * a4;
		Byeff = B.yval[i][j] * a1 + B.yval[i][j_plus] * a2 + B.yval[i + 1][j] * a3 + B.yval[i + 1][j_plus] * a4;
		Bzeff = B.zval[i][j] * a1 + B.zval[i][j_plus] * a2 + B.zval[i + 1][j] * a3 + B.zval[i + 1][j_plus] * a4;
		
		/* advance particle using modified Boris method */
		double beta = 1.0 / (1.0 + alpha2 * (Bxeff * Bxeff + Byeff * Byeff + Bzeff * Bzeff));
		//modifiedBorisPusher(species.part[n], alpha, beta, dt, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, loop);
		borisPusher(species.part[n], alpha, beta, dt, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff, loop);

		/* absorbing wall in x direction, periodic in y direction */
		if (species.part[i].x < lbound || species.part[i].x > rbound)
		{
			species.remParticle(i);
		}
		else if (species.part[i].y < dbound)
		{
			species.part[i].y += Ly;
		}
		else if (species.part[i].y > ubound)
		{
			species.part[i].y -= Ly;
		}
	}
}

void leapFrogRewind(Species& species, vectorField& E, vectorField& B, double rewindTime)
{
	/* effective electromagnetic field */
	double Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff;

	/* boris pusher coefficient */
	double alpha = 0.5 * species.specific_charge * rewindTime * lambda;  // lambda is the normalization factor
	double alpha2 = pow(alpha, 2);

	for (size_t n = 0; n < species.part.size(); n++)
	{
		double ratio_x = (species.part[n].x - lbound) * idx;
		int i = (int)floor(ratio_x);
		double prop_x1 = ratio_x - i;
		double prop_x2 = 1.0 - prop_x1;

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

		/* linear interploration */
		Exeff = E.xval[i][j] * a1 + E.xval[i][j_plus] * a2 + E.xval[i + 1][j] * a3 + E.xval[i + 1][j_plus] * a4;
		Eyeff = E.yval[i][j] * a1 + E.yval[i][j_plus] * a2 + E.yval[i + 1][j] * a3 + E.yval[i + 1][j_plus] * a4;
		Ezeff = E.zval[i][j] * a1 + E.zval[i][j_plus] * a2 + E.zval[i + 1][j] * a3 + E.zval[i + 1][j_plus] * a4;
		Bxeff = B.xval[i][j] * a1 + B.xval[i][j_plus] * a2 + B.xval[i + 1][j] * a3 + B.xval[i + 1][j_plus] * a4;
		Byeff = B.yval[i][j] * a1 + B.yval[i][j_plus] * a2 + B.yval[i + 1][j] * a3 + B.yval[i + 1][j_plus] * a4;
		Bzeff = B.zval[i][j] * a1 + B.zval[i][j_plus] * a2 + B.zval[i + 1][j] * a3 + B.zval[i + 1][j_plus] * a4;

		/* advance particle using modified Boris method */
		double beta = 1.0 / (1.0 + alpha2 * (Bxeff * Bxeff + Byeff * Byeff + Bzeff * Bzeff));
		modifiedBorisPusherRewind(species.part[n], alpha, beta, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff);
		borisPusherRewind(species.part[n], alpha, beta, Exeff, Eyeff, Ezeff, Bxeff, Byeff, Bzeff);

		/* absorbing wall in x direction, periodic in y direction */
		if (species.part[i].x < lbound || species.part[i].x > rbound)
		{
			species.remParticle(i);
		}
		else if (species.part[i].y < dbound)
		{
			species.part[i].y += Ly;
		}
		else if (species.part[i].y > ubound)
		{
			species.part[i].y -= Ly;
		}
	}
}