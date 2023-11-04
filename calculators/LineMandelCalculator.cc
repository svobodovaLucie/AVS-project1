/**
 * @file LineMandelCalculator.cc
 * @author Lucie Svobodová <xsvobo1x@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2023-10-27
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <immintrin.h>	// _mm_malloc()
#include <cstring>			// memcpy()

#include "LineMandelCalculator.h"

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int), 64));
	complex_tmp = (float *)(_mm_malloc(2 * width * sizeof(float), 64));
	tmp_data = (int *)(_mm_malloc(width * sizeof(int), 64));
	x_init = (float *)(_mm_malloc(width * sizeof(float), 64));
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(complex_tmp);
	complex_tmp = NULL;
	_mm_free(tmp_data);
	tmp_data = NULL;
	_mm_free(x_init);
	x_init = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	float *pcomplexReal;
	float *pcomplexImag;
	int *ptmp = tmp_data;
	float *pxInit = x_init;
	float y, r2, i2;
	int sum;

	// height loop
	for (int i = 0; i < height/2; i++) {
		// current imaginary value
		y = y_start + i * dy;

		// set the initial imaginary values
		std::uninitialized_fill(complex_tmp, complex_tmp + 2*width, y);
		// set the initial real values
		pcomplexReal = complex_tmp;
		pxInit = x_init;
		#pragma omp simd simdlen(64)
		for (int j = 0; j < width; j++) {
			*pcomplexReal = x_start + j * dx;
			*(pxInit++) = *pcomplexReal;
			pcomplexReal+=2;
		}

		// set the initial tmp_data with limit
		std::uninitialized_fill(tmp_data, tmp_data + width, limit);
		
		// iterations
		for (int l = 0; l < limit; l++) {
			pcomplexReal = complex_tmp;
			pcomplexImag = complex_tmp+1;
			ptmp = tmp_data;
			pxInit = x_init;
			sum = 0;

			// width loop 
			#pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
			for (int j = 0; j < width; j++) {
				// calculate one iteration of mandelbrot
				r2 = *pcomplexReal * *pcomplexReal;
				i2 = *pcomplexImag * *pcomplexImag;
				// update complex tmp
				*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
				*pcomplexReal = r2 - i2 + *pxInit;
				// update pointers
				pcomplexImag+=2;
				pcomplexReal+=2;
				pxInit++;

				if (r2 + i2 > 4.0f) {
					// update the current value if not updated yet
					(*ptmp != limit)?:*ptmp = l;
					sum++;
				}
				ptmp++;
			}
			// end the iterations before reaching the limit
			if (sum == width) break;
		}
		// update data 
		memcpy(data + (i*width), tmp_data, width*sizeof(int));
		memcpy(data + (height-i-1)*width, tmp_data, width*sizeof(int));
	}

	return data;
}
