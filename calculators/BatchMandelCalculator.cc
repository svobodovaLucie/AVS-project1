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
	tmp = (float *)(malloc(height * width * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	// @TODO cleanup the memory
	free(data);
	data = NULL;
	// free(tmp);
	// tmp = NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	int *pdata = data;
	float *ptmp = tmp;

	// initialize tmp buffer with x and y 
	for (int i = 0; i < height/2; i++) {
		for (int j = 0; j < width; j++) {
			*(ptmp) = x_start + j * dx;
			*(ptmp+1) = y_start + i * dy;
			ptmp = ptmp + 2;
		}
	}

	ptmp = tmp;

	const int N = limit; 
	const int TILE = 2;

	int cmp;
	int red = 0;
	float y;
	float x;
	float r2;
	float i2;

	for (int l = 0; l < limit; ++l) {

		for (int i = 0; i < height/2; i++)
		{
			y = y_start + i * dy; // current imaginary value

			#pragma omp simd
			for (int j = 0; j < width; j++)
			{
				x = x_start + j * dx; // current real value
				r2 = *(ptmp) * *(ptmp);
				i2 = *(ptmp+1) * *(ptmp+1);

				cmp = (r2 + i2 > 4.0f) ? 0 : 1;

				*(ptmp+1) = 2.0f * *(ptmp) * *(ptmp+1) + y;
				*(ptmp) = r2 - i2 + x;
				
				*(pdata) += cmp;
				
				// inverse pixel in matrix
				data[((height-i-1)*width)+j] += cmp;
				// data[i*width+j] += cmp;
				
				// update pointers
				pdata++;
				ptmp = ptmp + 2;
			}
		}
		pdata = data;
		ptmp = tmp;
	}
	return data;
}

