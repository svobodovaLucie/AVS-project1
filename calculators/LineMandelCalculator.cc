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


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	// @TODO allocate & prefill memory
	data = (int *)(malloc(height * width * sizeof(int)));
}

LineMandelCalculator::~LineMandelCalculator() {
	// @TODO cleanup the memory
	free(data);
	data = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	int *pdata = data;

	int half_height = height/2;
	int height_index, width_index, value;
	float x, y, zReal, zImag, r2, i2;

	for (int i = 0; i < half_height; i++) {
		y = y_start + i * dy;
		int height_offset = (height-i-1)*width;

		#pragma omp simd
		for (int j = 0; j < width; j++) {
			x = x_start + j * dx;

			// mandelbrot
			value = limit;
			zReal = x;
			zImag = y;

			for (int iteration = 0; iteration < limit; iteration++) {
				r2 = zReal * zReal;
				i2 = zImag * zImag;

				if (r2 + i2 > 4.0f) {
					value = iteration;
					break;
				}

				zImag = 2.0f * zReal * zImag + y;
				zReal = r2 - i2 + x;
			}
			*(pdata++) = data[height_offset + j] = value;
		}
	}
	return data;
}

/*
	float x, y, r2, i2;
	float zRealA[width];
	float zImagA[width];
	bool  b;
	float re, im;

	for (int i = 0; i < height/2; ++i) {
		y = y_start + i * dy; // current imaginary value	

		
		for (int l = 0; l < limit; ++l) {
			pdata = &data[i*width];
			// compute iterations for one row
			b = true;

			#pragma omp simd reduction (&:b)
			for (int j = 0; j < width; ++j) {

				x = x_start + j * dx; // current real value
				
				// nastaveni vychozi hodnoty, ze iterations == max
				if (l == 0) {
					*pdata = limit;
					data[((height-i-1)*width)+j] = limit;
					zRealA[j] = x;
					zImagA[j] = y;
				}

				re = zRealA[j];
				im = zImagA[j];

				// r2 = zRealA[j] * zRealA[j];
				// i2 = zImagA[j] * zImagA[j];
				r2 = re*re;
				i2 = im*im;

				if (r2 + i2 > 4.0f) {
					b &= false;

					// nastavit aktualni hloubku pouze tehdy, pokud jiz nebyla nastavena drive
					if (*pdata == limit) {
					*pdata = l;
					data[((height-i-1)*width)+j] = l;
					}
				}

				// zImagA[j] = 2.0f * zRealA[j] * zImagA[j] + y;
				// zRealA[j] = r2 - i2 + x;
				zImagA[j] = 2.0f * re * im + y;
				zRealA[j] = r2 - i2 + x;
				
				pdata++;
			}

			if (b) break;
		}
	}
	return data;
}
*/

