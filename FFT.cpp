/*
*  2N-point real-valued DFT/IDFT using N-point complex-valued FFT/IFFT.
*  no length check, no exception mechanism, sinTable & cosTable is precomputed for efficiency
*  edited by Huaxiang Zhang 2019/3/26
*/

#include "FFT.h"

/* private function for CFFT */ 
static size_t reverseBits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1U);
	return result;
}

/* CFFT implementation using Cooley-Tukey radix-2 method */
void CFFT(vector<double> &real, vector<double> &imag)
{   
	/* assuming real and imag have the same length */
	size_t n = real.size();

    /* compute levels = floor(log2(n)) */
	int levels = 0;  
	for (size_t temp = n; temp > 1U; temp >>= 1)
	{
        levels++;
	}

	/* bit-reversed addressing permutation */
	for (size_t i = 0; i < n; i++) 
	{
		size_t j = reverseBits(i, levels);
		if (j > i) 
		{
			swap(real[i], real[j]);
			swap(imag[i], imag[j]);
		}
	}

	/* Cooley-Tukey radix-2 method */
	for (size_t size = 2; size <= n; size *= 2) 
	{
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) 
		{
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) 
			{
				size_t l = j + halfsize;
				double tpre = real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
	    
		/* prevent overflow in 'size *= 2' */
	    if (size == n)  
		    break;
	}
}

void CIFFT(vector<double> &real, vector<double> &imag)
{
	/* swap the real part and imaginary part */
	CFFT(imag, real);

	/* normalization for inverse FFT */
	double norm = real.size();
	for (size_t i = 0; i < norm; i++)
	{
		real[i] /=  norm;
		imag[i] /=  norm;
	}
}

void RDFT(vector<double> &inputReal, vector<double> &real, vector<double> &imag)
{
	/* divide original sequence into even and odd subsequence */
	int subLength = inputReal.size() / 2;
	vector<double> even(subLength );
	vector<double> odd(subLength );

	for (int i = 0; i < subLength; i++)
	{
		even[i] = inputReal[2 * i];
		odd[i] = inputReal[2 * i + 1];
	}
	
	/* using complex FFT */
	CFFT(even, odd);

	for (int i = 1; i < subLength; i++)
	{
		real[i] = even[i] * ArTable[i] - odd[i] * AiTable[i] + even[subLength - i] * BrTable[i] + odd[subLength - i] * BiTable[i];
		real[2 * subLength - i ] = real[i];

		imag[i] = odd[i] * ArTable[i] + even[i] * AiTable[i] + even[subLength - i] * BiTable[i] - odd[subLength - i] * BrTable[i];
		imag[2 * subLength - i ] = -1 * imag[i];
	}

	real[0] = even[0] + odd[0];
	imag[0] = 0;
	real[subLength] = even[0] - odd[0];
	imag[subLength] = 0;
}

void RIDFT(vector<double> &outputReal, vector<double> &real, vector<double> &imag)
{
	/* divide original sequence into even and odd subsequence */
	int subLength = outputReal.size() / 2;
	vector<double> even(subLength);
	vector<double> odd(subLength);

	for (int i = 1; i < subLength; i++)
	{
		even[i] = real[i] * ArTable[i] + imag[i] * AiTable[i] + real[subLength - i] * BrTable[i] - imag[subLength - i] * BiTable[i];

		odd[i] = imag[i] * ArTable[i] - real[i] * AiTable[i] - real[subLength - i] * BiTable[i] - imag[subLength - i] * BrTable[i];
	}

	even[0] = (real[0] + real[subLength]) / 2;
	odd[0] = (real[0] - real[subLength]) / 2;

	/* using complex FFT */
	CIFFT(even, odd);

	for (int i = 0; i < subLength; i++)
	{
		outputReal[2 * i] = even[i];
		outputReal[2 * i + 1] = odd[i];
	}
}