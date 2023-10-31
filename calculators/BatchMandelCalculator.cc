/**
 * @file BatchMandelCalculator.cc
 * @author Lucie Svobodová <xsvobo1x@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 2023-10-27
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	// @TODO allocate & prefill memory
	data = (int *)(malloc(height * width * sizeof(int)));
	// tmp = (float *)(malloc(2 * height * width * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	// @TODO cleanup the memory
	free(data);
	data = NULL;
	// free(tmp);
	// tmp = NULL;
}

/*
int * BatchMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	// @TODO implement the calculator & return array of integers
	int *pdata = data;

	float x, y, r2, i2;
	float zRealA[width];
	float zImagA[width];
	bool  b;
	float re, im;

	int N = 0;
	int TILE = 0;

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
					// tmp[j+j] = x;
					// tmp[j+j+1] = y;
				}

				re = zRealA[j];
				im = zImagA[j];
				// re = tmp[j+j];
				// im = tmp[j+j+1];

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

				zImagA[j] = 2.0f * zRealA[j] * zImagA[j] + y;
				zRealA[j] = r2 - i2 + x;
				// tmp[j+j+1] = 2.0f * re * im + y;
				// tmp[j+j] = r2 - i2 + x;
				
				pdata++;
			}

			if (b) break;
		}
	}
	return data;
*/

int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data;

	int size = height / 2;
	int TILE = 64;	// TODO
	int p = 2;	// TODO parametr, modulo, celociselne deleni atd.
	int N_width = width/TILE;
	int N_height = height/(TILE*2);
	bool b;

	// height
	for (int i = 0; i < N_height; i++) {
		for (int ii = 0; ii < TILE; ii++) {

			// width
			for (int j = 0; j < N_width; j++) {
				#pragma omp simd
				for (int jj = 0; jj < TILE; jj++) {
					int width_index = j * TILE + jj;
					int height_index = i * TILE + ii;

					float x = x_start + width_index * dx;
					float y = y_start + height_index * dy;

					// mandelbrot
					int value = limit;
					float zReal = x;
					float zImag = y;

					// b = true;
					// #pragma omp simd // reduction (&:b)
					for (int iteration = 0; iteration < limit; iteration++) {
						float r2 = zReal * zReal;
						float i2 = zImag * zImag;

						if (r2 + i2 > 4.0f) {
							value = iteration;
							break;
							// b &= false;
							// b = false;
							// if (value == limit) {
							// 	value = iteration;
							// }
						}

						zImag = 2.0f * zReal * zImag + y;
						zReal = r2 - i2 + x;
					}
					pdata[height_index * width + width_index] = value;
					pdata[(height - height_index - 1) * width + width_index] = value;
					// if (b) break;
				}
			}
		}
	}
	return data;
}

