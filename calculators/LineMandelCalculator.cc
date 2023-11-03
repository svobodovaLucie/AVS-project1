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
#include <cstring>			// memset()

#include "LineMandelCalculator.h"

#define ALIGN 64

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	//data = (int *)(malloc(height * width * sizeof(int)));
	data = (int *)(_mm_malloc(height * width * sizeof(int), ALIGN));
	complex_tmp = (float *)(_mm_malloc(2 * width * sizeof(float), ALIGN));
	x_real = (float *)(_mm_malloc(width * sizeof(float), ALIGN));
	tmp_data = (int *)(_mm_malloc(width * sizeof(int), ALIGN));
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(complex_tmp);
	complex_tmp = NULL;
	_mm_free(x_real);
	x_real = NULL;
	_mm_free(tmp_data);
	complex_tmp = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *pcomplexReal = complex_tmp;
	float *pcomplexImag = complex_tmp+1;
	float *pxReal = x_real;
	int *ptmp = tmp_data;
	float y, x, r2, i2;

	// initialize the complex_tmp array
	for (int i = 0; i < height/2; i++) {

		y = y_start + i * dy;				// current imaginary value

		// set the initial values
		std::uninitialized_fill(complex_tmp, complex_tmp + 2*width, y);
		// set the initial values for x
		pcomplexReal = complex_tmp;
		pcomplexImag = complex_tmp+1;
		pxReal = x_real;

		#pragma omp simd simdlen(64)
		for (int j = 0; j < width; j++) {
			*pxReal = x_start + j * dx;
			*pcomplexReal = *(pxReal++);
			pcomplexReal+=2;
		}

		std::uninitialized_fill(tmp_data, tmp_data + width, limit);
		ptmp = tmp_data;
		
		// iterations
		for (int l = 0; l < limit; l++) {
			pcomplexReal = complex_tmp;
			pcomplexImag = complex_tmp+1;
			pxReal = x_real;
			ptmp = tmp_data;

			int sum = 0;
			#pragma omp simd reduction(+:sum) simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) // simdlen(1024)
			for (int j = 0; j < width; j++) {

				//x = x_start + j * dx;		// current real value

				// calculate one iteration of mandelbrot
				// with the use of tmp value in pcomplex
				r2 = *pcomplexReal * *pcomplexReal;
				i2 = *pcomplexImag * *pcomplexImag;

				if (r2 + i2 > 4.0f) {
					if (*ptmp == limit) {
						*ptmp = l;
					}
					sum++;
				}

				// update complex tmp
				*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
				*pcomplexReal = r2 - i2 + *(pxReal++);
	
				pcomplexImag+=2;
				pcomplexReal+=2;

				ptmp++;
			}

			// if all are good -> end the iterations loop
			if (sum == width) break;
		}
		memcpy(data + (i*width), tmp_data, width*sizeof(int));
		memcpy(data + (height-i-1)*width, tmp_data, width*sizeof(int));
	}

	return data;
}
