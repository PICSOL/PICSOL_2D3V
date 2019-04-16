#include "Global.h"

/*
Global variables definition
*/

Param param;
const size_t rank_x = param.GRID_NUM_X;
const size_t rank_y = param.GRID_NUM_Y;
const double lbound = param.LBOUND;
const double rbound = param.RBOUND;
const double ubound = param.UBOUND;
const double dbound = param.DBOUND;
const double dx = param.DELTA_X;
const double idx = 1.0 / dx;
const double dx2 = pow(param.DELTA_X, 2);
const double idx2 = 1.0 / dx2;
const double dy = param.DELTA_Y;
const double idy = 1.0 / dy;
const double dy2 = pow(param.DELTA_Y, 2);
const double idy2 = 1.0 / dy2;
const double dxy2 = -2.0 * (1.0 / dy2 + 1.0 / dx2);
const double idxy2 = 1.0 / dxy2;
const double Lx = dx * (rank_x - 1);
const double Ly = dy * rank_y;
const double twoPi = 2.0 * param.PI;
const size_t maxIteration = param.MAX_ITERATION;
const double minError = param.MIN_ERROR;
const double ioncharge = param.IONCHARGE;
const double electcharge = param.ELECTCHARGE;
const int initParticle = param.INITIAL_PARTICLE;

const double lambda = param.EPS_0 * param.MU_0 * param.ETEMP * param.QE / param.ME;
const double normEl = sqrt(param.ETEMP * param.QE * param.DENSITY / param.EPS_0);
const double normEt = normEl * lambda;
const double normB = normEl * param.EPS_0 * param.MU_0 * sqrt(param.ETEMP * param.QE / param.ME);

const Vector3D<double> B_0(param.BX/normB, param.BY/normB, 0);
vector<double> omegaTable(2 * param.MAX_ITERATION);
vector<double> phi_m_real(param.GRID_NUM_Y);
vector<double> phi_m_imag(param.GRID_NUM_Y);
vector<double> psi_m_real(param.GRID_NUM_Y);
vector<double> psi_m_imag(param.GRID_NUM_Y);
vector<double> FFTbound(param.GRID_NUM_Y);
vector<double> diagElement(param.GRID_NUM_Y / 2 + 1);
vector<double> ArTable(param.GRID_NUM_Y / 2);
vector<double> AiTable(param.GRID_NUM_Y / 2);
vector<double> BrTable(param.GRID_NUM_Y / 2);
vector<double> BiTable(param.GRID_NUM_Y / 2);
vector<double> sinTable(param.GRID_NUM_Y / 4);
vector<double> cosTable(param.GRID_NUM_Y / 4);

void initialGlobal()
{
	double Pi = param.PI;
	double rho_Jac = (cos(Pi / rank_x) + pow(dx / dy, 2)*cos(2 * Pi / rank_y)) / (1 + pow(dx / dy, 2));
	
	omegaTable[0] = 1;
	omegaTable[1] = 1.0 / (1 - 0.5 * pow(rho_Jac, 2));

	for (int i = 2; i < 2*param.MAX_ITERATION; i++)
	{
		omegaTable[i] = 1.0 / (1 - pow(rho_Jac, 2)*omegaTable[i - 1] / 4.0);
	}

	for (int i = 0; i <= param.GRID_NUM_Y / 2; i++)
	{
		diagElement[i] = -1 * (2 + pow(twoPi*i / Ly, 2)*pow(dx, 2));
	}

	for (int i = 0; i < param.GRID_NUM_Y / 4; i++)
	{
		cosTable[i] = cos(4 * Pi * i / param.GRID_NUM_Y);
		sinTable[i] = sin(4 * Pi * i / param.GRID_NUM_Y);
	}

	for (int i = 0; i < param.GRID_NUM_Y / 2; i++)
	{
		double cosk = cos(2 * Pi * i / param.GRID_NUM_Y);
		double sink = sin(2 * Pi * i / param.GRID_NUM_Y);
		ArTable[i] = 0.5 * (1.0 - sink);
		AiTable[i] = -0.5 * cosk;
		BrTable[i] = 0.5 * (1.0 + sink);
		BiTable[i] = 0.5 * cosk;
	}

	RDFT(param.LB, phi_m_real, phi_m_imag);
	RDFT(param.RB, psi_m_real, psi_m_imag);
}