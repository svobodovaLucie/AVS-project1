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
#include <immintrin.h>
#include <cstring>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int), 64));
	complex_tmp = (float *)(_mm_malloc(2 * 64 * sizeof(float), 64));
	// tmp_data = (int *)(_mm_malloc(64 * sizeof(int), 64));
	x_init = (float *)(_mm_malloc(64 * sizeof(float), 64));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(complex_tmp);
	complex_tmp = NULL;
	// _mm_free(tmp_data);
	// complex_tmp = NULL;
	_mm_free(x_init);
	x_init = NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {
	float *pcomplexReal;
	float *pcomplexImag;
	float *pxInit = x_init;
	float y, x, r2, i2;
	int sum;
	int TILE = 64;
	int i_width, height_i_1, j_TILE;

	int N_width = width/TILE;
	int TILE_sizeofint = TILE*sizeof(int);

	std::uninitialized_fill(data, data + width*height, limit);

	// height loop
	for (int i = 0; i < height/2; i++) {

		// current imaginary value
		y = y_start + i * dy;

		i_width = i*width;
		height_i_1 = (height-i-1)*width;

		#pragma omp simd simdlen(64)
		for (int j = 0; j < N_width; j++) {
				
			// set the initial imaginary values
			std::uninitialized_fill(complex_tmp, complex_tmp + 2*TILE, y);
			// set the initial real values
			pcomplexReal = complex_tmp;
			pxInit = x_init;
			j_TILE = j*TILE;

			#pragma omp simd simdlen(64) aligned(pcomplexReal,pxInit:64)
			for (int jj = 0; jj < TILE; jj++) {
				*pcomplexReal = x_start + (j_TILE+jj) * dx;
				*(pxInit++) = *pcomplexReal;
				pcomplexReal+=2;
			}

			
			// iterations
			for (int l = 0; l < limit; l++) {
				pcomplexReal = complex_tmp;
				pcomplexImag = complex_tmp+1;
				pxInit = x_init;
				sum = 0;

				// width loop 
				#pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
				// inner loop with TILE size
				for (int jj = 0; jj < TILE; jj++) {
					// calculate one iteration of mandelbrot
					r2 = *pcomplexReal * *pcomplexReal;
					i2 = *pcomplexImag * *pcomplexImag;

					if (r2 + i2 > 4.0f) {
						if (data[i_width+j_TILE+jj] == limit) {
							data[i_width+j_TILE+jj] = l;
							data[height_i_1+j_TILE+jj] = l;
						}
						sum++;
					}

					// update complex tmp
					*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
					*pcomplexReal = r2 - i2 + *pxInit;
		
					pcomplexImag+=2;
					pcomplexReal+=2;
					pxInit++;
				}
				// end the iterations before reaching the limit
				if (sum == width) break;
			}
		}
	}
	return data;
}