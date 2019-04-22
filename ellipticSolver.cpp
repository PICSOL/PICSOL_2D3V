#include "ellipticSolver.h"

void tridiagSolver(double b, vector<double> &x, vector<double> &y)
{   
	int n = x.size();
	vector<double> c_prime(n - 1);

	/* initialization, storing d_prime into x */
	c_prime[0] = 1.0 / b;
	x[0] = y[0] / b;

	/* forward subsititution */
	for (int i = 1; i < n - 1; i++)
	{
		c_prime[i] = 1.0 / (b - c_prime[i - 1]);
		x[i] = c_prime[i] * (y[i] - x[i - 1]);
	}
	x[n - 1] = (y[n - 1] - x[n - 2]) / (b - c_prime[n - 2]);

	/* backward subsititution */
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] -= c_prime[i] * x[i + 1];
	}
}

/* in this function index 'i' stands for x direction, index 'j' stands for y direction */
void poissonFFT(vector<vector<double>> &P, vector<vector<double>> &Q, vector<double> &LBm_real,
                vector<double> &LBm_imag, vector<double> &RBm_real, vector<double> &RBm_imag)
{
   /* real and imaginary part are computed separatedly */
	vector<vector<double>> Qm_real, Qm_imag, Pm_real, Pm_imag;
    vector<double> real(rank_y), imag(rank_y);
	vector<double> qx_real(rank_x - 2),qx_imag(rank_x - 2);
	vector<double> px_real(rank_x - 2),px_imag(rank_x - 2);

	/* loop for calling FFT in y direction for Qm */
	for (size_t i = 1; i < rank_x-1; i++)
	{
		/* real-valued FFT */
		RDFT(Q[i], real, imag);
		Qm_real.push_back(real);
		Qm_imag.push_back(imag);
	}

	/* Physical quantities are real numbers, only need to compute from 0 to rank_y / 2 */
	for (size_t j = 0; j <= rank_y / 2; j++)
	{
		/* target vector formation */
		for (size_t i = 0; i < rank_x - 2; i++)
		{
			qx_real[i] = Qm_real[i][j] * dx2;
			qx_imag[i] = Qm_imag[i][j] * dx2;
		}

		/* impose dirichlet boundary condition */
		qx_real[0] -= LBm_real[j];
		qx_imag[0] -= LBm_imag[j];
		qx_real[rank_x - 3] -= RBm_real[j];
		qx_imag[rank_x - 3] -= RBm_imag[j];

		/* solve tridiagnol matrix */
		tridiagSolver(diagElement[j], px_real, qx_real);
		tridiagSolver(diagElement[j], px_imag, qx_imag);

		/* store the Pm matrix */
		Pm_real.push_back(px_real);
		Pm_imag.push_back(px_imag);
	}

	/* loop for calling IFFT in y direction for P */
	for (size_t i = 1; i < rank_x-1; i++)
	{
		/* target vector formation */
		for (size_t j = 0; j <= rank_y / 2; j++)
		{
			real[j] = Pm_real[j][i - 1];
			imag[j] = Pm_imag[j][i - 1];
		}

		/* real-valued inverse FFT */
		RIDFT(P[i], real, imag);
	}
}

void poissonSOR(vector<vector<double>> &P, vector<vector<double>> &Q)
{   
	/* iteration criterion */
	double Qnorm = 0;
	for (size_t i = 1; i < rank_x - 1; i++) 
	{
		for (size_t j = 0; j < rank_y; j++) 
		{
             Qnorm += abs(Q[i][j]);
        }
	}		
	Qnorm *= minError;

	for (size_t iter = 0; iter < maxIteration; iter++) 
	{
		/* the sum of residual on every grid point */
		double norm = 0;
		double residual;

		/* chebyshev acceleration SOR coefficient */
		double omega = omegaTable[iter];

		/* dirichlet boundary condition in x direction */
		for (size_t i = 1; i < rank_x - 1; i++)
		{
			residual = idx2 * (P[i + 1][0] + P[i - 1][0])
				+ idy2 * (P[i][1] + P[i][rank_y - 1])
				+ dxy2 * P[i][0]
				- Q[i][0];

			norm += abs(residual);
			P[i][0] -= omega * residual * idxy2;
			
			/* periodic boundary condition in y direction */
			for (size_t j = 1; j < rank_y - 1; j ++)
			{
				residual = idx2 * (P[i + 1][j] + P[i - 1][j])
					     + idy2 * (P[i][j+1] + P[i][j-1])
					     + dxy2 * P[i][j]
					     - Q[i][j];

				norm += abs(residual);
				P[i][j] -= omega * residual * idxy2;
			}

			residual = idx2 * (P[i + 1][rank_y - 1] + P[i - 1][rank_y - 1])
				+ idy2 * (P[i][rank_y - 2] + P[i][0])
				+ dxy2 * P[i][rank_y - 1]
				- Q[i][rank_y - 1];

			norm += abs(residual);
			P[i][rank_y - 1] -= omega * residual * idxy2;
		}
		if (norm < Qnorm)
			return;
	}
	cerr << "poissonSOR reaches maximum iterations !" << endl;
}

void helmholtzSOR(vector<vector<double>> &P, vector<vector<double>> &Q, vector<vector<double>> &a)
{
	/* iteration criterion */
	double Qnorm = 0;
	for (size_t i = 1; i < rank_x - 1; i++)
	{
		for (size_t j = 0; j < rank_y; j++)
		{
			Qnorm += abs(Q[i][j]);
		}
	}
	Qnorm *= minError;

	for (size_t iter = 0; iter < maxIteration; iter++)
	{
		/* the sum of residual on every grid point */
		double norm = 0;
		double residual;

		/* chebyshev acceleration SOR coefficient */
		double omega = omegaTable[iter];

		/* dirichlet boundary condition in x direction */
		for (size_t i = 1; i < rank_x - 1; i++)
		{
			double e = dxy2 + a[i][0];
			residual = idx2 * (P[i + 1][0] + P[i - 1][0])
				+ idy2 * (P[i][1] + P[i][rank_y - 1])
				+ e * P[i][0]
				- Q[i][0];

			norm += abs(residual);
			P[i][0] -= omega * residual / e;

			/* periodic boundary condition in y direction */
			for (size_t j = 1; j < rank_y - 1; j++)
			{
				e = dxy2 + a[i][j];
				residual = idx2 * (P[i + 1][j] + P[i - 1][j])
					+ idy2 * (P[i][j + 1] + P[i][j - 1])
					+ e * P[i][j]
					- Q[i][j];

				norm += abs(residual);
				P[i][j] -= omega * residual / e;
			}

			e = dxy2 + a[i][rank_y - 1];
			residual = idx2 * (P[i + 1][rank_y - 1] + P[i - 1][rank_y - 1])
				+ idy2 * (P[i][rank_y - 2] + P[i][0])
				+ e * P[i][rank_y - 1]
				- Q[i][rank_y - 1];

			norm += abs(residual);
			P[i][rank_y - 1] -= omega * residual / e;
		}
		if (norm < Qnorm)
			return;
	}
	cerr << "helmholtzSOR reaches maximum iterations !" << endl;
}