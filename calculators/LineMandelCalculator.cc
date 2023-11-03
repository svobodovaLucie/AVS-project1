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
	complexReal = (float *)(_mm_malloc(width * sizeof(float), ALIGN));
	complexImag = (float *)(_mm_malloc(width * sizeof(float), ALIGN));
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(complex_tmp);
	complex_tmp = NULL;
	_mm_free(complexReal);
	complex_tmp = NULL;
	_mm_free(complexImag);
	complex_tmp = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	// float *pcomplexReal = complex_tmp;
	// float *pcomplexImag = complex_tmp+1;
	float *real = complexReal;
	float *imag = complexImag;

	float y, x, r2, i2;

	// initialize data array to limit value
	std::uninitialized_fill(data, data + width*height, limit);

	// initialize the complex_tmp array
	for (int i = 0; i < height/2; i++) {
		y = y_start + i * dy;				// current imaginary value

		// set the initial values
		// std::uninitialized_fill(complex_tmp, complex_tmp + 2*width, y);
		// std::uninitialized_fill(complexReal, complex_tmp + width, x);
		std::uninitialized_fill(complexImag, complexImag + width, y);
		// set the initial values for x

		real = complexReal;
		#pragma omp simd aligned(real: 64)
		for (int j = 0; j < width; j++) {
			*(real++) = x_start + j * dx;
			// *pcomplexReal = x;
			// pcomplexReal += 2;
		}

		// pcomplexReal = complex_tmp;
		// pcomplexImag = complex_tmp+1;

		// pcomplex = complex_tmp;
		// iterations
		for (int l = 0; l < limit; l++) {

			// width loop
			//bool b = false;
			real = complexReal;
			imag = complexImag;

			int sum = 0;
			#pragma omp simd reduction(+:sum) aligned(real, imag: 64) // simdlen(1024)
			for (int j = 0; j < width; j++) {
				// stop the width loop when reduction...

				x = x_start + j * dx;		// current real value

				// calculate one iteration of mandelbrot
				// with the use of tmp value in pcomplex
				// TODO transform to pointer++
				// zReal = complex_tmp[2*j];
				// zImag = complex_tmp[2*j+1];

				// zReal = *pcomplexReal;
				// zImag = *pcomplexImag;

				// r2 = zReal * zReal;
				// i2 = zImag * zImag;
				// r2 = *pcomplexReal * *pcomplexReal;
				// i2 = *pcomplexImag * *pcomplexImag;
				r2 = *real * *real;
				i2 = *imag * *imag;

				if (r2 + i2 > 4.0) {
					if (data[i*width+j] == limit) {
					//if (*pdata == -1) {
						//*pdata = l;
						data[i*width+j] = l;
						// *pdata = l;
						data[(height-i-1)*width+j] = l;
					}
					sum++;
				}

				// update complex tmp
				// complex_tmp[2*j+1] = 2.0f * zReal * zImag + y;
				// complex_tmp[2*j] = r2 - i2 + x;
				// *pcomplexImag = 2.0f * zReal * zImag + y;
				// *pcomplexReal = r2 - i2 + x;
				// *pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
				// *pcomplexReal = r2 - i2 + x;
				*imag = 2.0f * *real * *imag + y;
				*real = r2 - i2 + x;

				// pcomplexImag+=2;
				// pcomplexReal+=2;
				real++;
				imag++;
			}
			// pcomplexReal = complex_tmp;
			// pcomplexImag = complex_tmp+1;

			// if all are good -> end the iterations loop
			if (sum == width) break;
		}
		// pcomplexImag = complex_tmp+1;
		// pcomplexReal = complex_tmp;
		// update pdata pointer
	}

	// for (int i = 0; i < width*height; ++i) {
	// 	printf("%d ", data[i]);
	// }

	real = complexReal;
	imag = complexImag;

	return data;
}
