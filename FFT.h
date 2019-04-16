/*
*  2N-point real-valued DFT/IDFT using N-point complex-valued FFT/IFFT.
*  no length check, no exception mechanism, sinTable & cosTable is precomputed for efficiency
*  edited by Huaxiang Zhang 2019/3/26
*/

#pragma once
#include "Global.h"
using namespace std;

/*
Coefficients are defined as below
For CFFT/ CIFFT

    for (int i = 0; i < n / 2; i++) 
    {
        cosTable[i] = cos( 2 * Pi * i / n);
        sinTable[i] = sin( 2 * Pi * i / n);
    }

For RDFT/RIDFT

	for (int i = 0; i < n / 4; i++) 
	{
		cosTable[i] = cos( 4 * Pi * i / n);
		sinTable[i] = sin( 4 * Pi * i / n);
	}
	
	for (int i = 0; i < n / 2; i++) 
	{
		double cosk = cos(2 * Pi * i / n);
		double sink = sin(2 * Pi * i / n);
		ArTable[i] = 0.5 * (1.0 - sink);
		AiTable[i] = -0.5 * cosk;
		BrTable[i] = 0.5 * (1.0 + sink);
		BiTable[i] = 0.5 * cosk;
	}
*/

/*
Computes complex-valued FFT of the given complex vector, storing the result back into the vector container. 
The vector must have length which is the power of 2, using the Cooley-Tukey radix-2 algorithm.
Coefficient sinTable, cosTable is precomputed 
*/
void CFFT(vector<double> &real, vector<double> &imag);

/*
Computes complex-valued inverse FFT of the given complex vector, storing the result back into the vector container. 
The vector must have length which is the power of 2, using the Cooley-Tukey radix-2 algorithm.
Coefficient sinTable, cosTable is precomputed 
*/
void CIFFT(vector<double> &real, vector<double> &imag);

/*
Computes real-valued DFT of given real vector, storing the result back into the vector.
The vector must have length which is the power of 2, using N-point complex-valued FFT to implement 2N-point real-valued DFT.
Coefficient sinTable, cosTable, Ar, Ai, Br, Bi is precomputed 
*/
void RDFT(vector<double> &inputReal, vector<double> &real, vector<double> &imag);

/*
Computes real-valued inverse DFT of given real vector, storing the result back into the vector.
The vector must have length which is the power of 2, using N-point complex-valued inverse FFT to implement 2N-point real-valued inverse DFT.
coefficient sinTable, cosTable, Ar, Ai, Br, Bi is precomputed 
*/
void RIDFT(vector<double> &outputReal, vector<double> &real, vector<double> &imag);



