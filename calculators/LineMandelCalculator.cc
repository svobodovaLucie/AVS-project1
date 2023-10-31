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

	float x, y, r2, i2;
	float zRealA[width];
	float zImagA[width];
	bool  c[width];
	bool  b;
	int l_history;

	for (int i = 0; i < height/2; ++i) {
		y = y_start + i * dy; // current imaginary value		
		for (int l = 0; l < limit; ++l) {
			pdata = &data[i*width];
			// compute iterations for one row
			// b = (l == 0) ? true : b;
b = true;
			// #pragma omp simd
			#pragma omp simd reduction (&:b)
			for (int j = 0; j < width; ++j) {
				// c uz ma dopocitanou max iteraci
				// c[j] = (l == 0) ? false : c[j];
				
				// nastaveni vychozi hodnoty, ze iterations == max
				if (l == 0) {
					*pdata = limit;
					data[((height-i-1)*width)+j] = limit;
				}
				
				// TODO b
				// if (c[j]) {
				// 	pdata++;
				// 	continue;
				// }

				// one pixel
				x = x_start + j * dx; // current real value

				zRealA[j] = (l == 0) ? x : zRealA[j];
				zImagA[j] = (l == 0) ? y : zImagA[j];

				r2 = zRealA[j] * zRealA[j];
				i2 = zImagA[j] * zImagA[j];

				l_history = *pdata;

				if (r2 + i2 > 4.0f) {
					b &= false;
					// nastavit aktualni hloubku pouze tehdy, pokud jiz nebyla nastavena drive
					if (*pdata == limit) {
					*pdata = l;
					data[((height-i-1)*width)+j] = l;
					}
				}

				pdata++;

				zImagA[j] = 2.0f * zRealA[j] * zImagA[j] + y;
				zRealA[j] = r2 - i2 + x;
				
			}
			if (b) {
				break;
			}
		}
	}
	return data;
}
