#include "fieldUpdater.h"

/* temporary storage for RHS(right hand side of equaitons) */
static scalarField sTemp(param);
static vectorField vTemp(param);
static scalarField Theta(param);


void rhoComp(scalarField& electDen, scalarField& ionDen, scalarField& rho)
{
	for (size_t i = 0; i < rank_x; i++)
	{
		for (size_t j = 0; j < rank_y; j++)
		{
             rho.val[i][j] = ioncharge * ionDen.val[i][j] + electcharge * electDen.val[i][j];
		}		
	}
}

void JComp(vectorField& electCurrent, vectorField& ionCurrent, vectorField& J)
{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	J.sum(electCurrent, ionCurrent);
}

void convectComp(vectorField& electCon, vectorField& ionCon, vectorField& convect)
{
	convect.sum(electCon, ionCon);
}

void EtotalComp(vectorField& E_l, vectorField& E_t, vectorField& E_total)
{
	for (size_t i = 0; i < rank_x; i++)
	{
		for (size_t j = 0; j < rank_y; j++)
		{
			E_total.xval[i][j] = E_t.xval[i][j] + E_l.xval[i][j] / lambda;
			E_total.yval[i][j] = E_t.yval[i][j] + E_l.yval[i][j] / lambda;
			E_total.zval[i][j] = E_t.zval[i][j];
		}
	}
}

void phiComp(scalarField &rho, scalarField &phi)
{
	for (size_t i = 0; i < rank_x; i++)
	{
		for (size_t j = 0; j < rank_y; j++)
		{
			sTemp.val[i][j] = - rho.val[i][j];
		}
	}
	
	//poissonSOR(phi.val, sTemp.val);
	poissonFFT(phi.val, sTemp.val, phi_m_real, phi_m_imag, psi_m_real, psi_m_imag);
}

void ElComp(scalarField &phi, vectorField &E_l)
{
	
	gradient(E_l, phi, -1.0);

	/* forward and backward difference at conducting wall boundaries */
	for (size_t j = 1; j < rank_y - 1; j++)
	{
		E_l.xval[0][j] = idx * (phi.val[0][j] - phi.val[1][j]);
		E_l.xval[rank_x - 1][j] = idx * (phi.val[rank_x - 2][j] - phi.val[rank_x - 1][j]);
	    
		// y-direction is questionable
		E_l.yval[0][j] = 0.5 * idy * (phi.val[0][j - 1] - phi.val[0][j + 1]);
		E_l.yval[rank_x - 1][j] = 0.5 * idy * (phi.val[rank_x - 1][j - 1] - phi.val[rank_x - 1][j + 1]);
	}

	/* periodic in y direction */
	E_l.xval[0][0] = idx * (phi.val[0][0] - phi.val[1][0]);
	E_l.xval[rank_x - 1][0] = idx * (phi.val[rank_x - 2][0] - phi.val[rank_x - 1][0]);
	E_l.xval[0][rank_y - 1] = idx * (phi.val[0][rank_y - 1] - phi.val[1][rank_y - 1]);
	E_l.xval[rank_x - 1][rank_y - 1] = idx * (phi.val[rank_x - 2][rank_y - 1] - phi.val[rank_x - 1][rank_y - 1]);

	// y-direction is questionable
	E_l.yval[0][0] = 0.5 * idy * (phi.val[0][rank_y - 1] - phi.val[0][1]);
	E_l.yval[rank_x - 1][0] = 0.5 * idy * (phi.val[rank_x - 1][rank_y - 1] - phi.val[rank_x - 1][1]);
	E_l.yval[0][rank_y - 1] = 0.5 * idy * (phi.val[0][rank_y - 2] - phi.val[0][0]);
	E_l.yval[rank_x - 1][rank_y - 1] = 0.5 * idy * (phi.val[rank_x - 1][rank_y - 2] - phi.val[rank_x - 1][0]);
} 

void BComp(vectorField &J, vectorField &B)
{
	rotation(vTemp, J, -1.0);

	/* OpenMP multithreading area for field updater */
    #pragma omp parallel default(shared)
	{
        #pragma omp sections nowait
		{
            #pragma omp section
			poissonFFT(B.xval, vTemp.xval, FFTbound, FFTbound, FFTbound, FFTbound);
			//poissonSOR(B.xval, vTemp.xval);

            #pragma omp section
			poissonFFT(B.yval, vTemp.yval, FFTbound, FFTbound, FFTbound, FFTbound);
			//poissonSOR(B.yval, vTemp.yval);

            #pragma omp section
			poissonFFT(B.zval, vTemp.zval, FFTbound, FFTbound, FFTbound, FFTbound);
			//poissonSOR(B.zval, vTemp.zval);
		}
	}

	/* add background field */
	B.add(B_0);
}

void EtComp(scalarField &rho, vectorField &Current, vectorField &Convect, vectorField &E_l, vectorField &B, vectorField &E_t)
{
	/* pre-Computaion process for A (coefficient matrix) and Q(right hand side) in helmholtz equation */
	for (size_t i = 0; i < rank_x ; i++)
	{
		for (size_t j = 0; j < rank_y ; j++)
		{
			sTemp.val[i][j] = - rho.val[i][j] * lambda;
			vTemp.xval[i][j] = Convect.xval[i][j] + rho.val[i][j] * E_l.xval[i][j]
				         + lambda * (Current.yval[i][j] * B.zval[i][j] - Current.zval[i][j] * B.yval[i][j]);
			vTemp.yval[i][j] = Convect.yval[i][j] + rho.val[i][j] * E_l.yval[i][j]
				         + lambda * (Current.zval[i][j] * B.xval[i][j] - Current.xval[i][j] * B.zval[i][j]);
			vTemp.zval[i][j] = Convect.zval[i][j] 
				         + lambda * (Current.xval[i][j] * B.yval[i][j] - Current.yval[i][j] * B.xval[i][j]);
		}
	}

	/* OpenMP multithreading area for field updater */
    #pragma omp parallel default(shared)
	{
        #pragma omp sections nowait
		{
            #pragma omp section
			helmholtzSOR(E_t.xval, vTemp.xval, sTemp.val);

            #pragma omp section
			helmholtzSOR(E_t.yval, vTemp.yval, sTemp.val);

            #pragma omp section
			helmholtzSOR(E_t.zval, vTemp.zval, sTemp.val);
		}
	}

	divergence(sTemp, E_t, 1.0);

	poissonFFT(Theta.val, sTemp.val, FFTbound, FFTbound, FFTbound, FFTbound);	
	//poissonSOR(Theta.val, sTemp.val);
    
	/* substract the divergence part to abtain authentic E_t */
	gradient(vTemp, Theta, 1.0);
	E_t.substract(vTemp);
}

