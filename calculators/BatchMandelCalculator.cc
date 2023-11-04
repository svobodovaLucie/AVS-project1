/**
 * @file BatchMandelCalculator.cc
 * @author Lucie Svobodová <xsvobo1x@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 2023-11-03
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdexcept>
#include <immintrin.h>

#include "BatchMandelCalculator.h"

#define TILE 64

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int), TILE));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	data = NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {

	float x,y,r2,i2,zReal,zImag;
	int height_index, width_index, value;
	int *pdata = data;
	int width_tile = width/TILE;
	int tiles_half = (width*height)/(TILE*TILE*2);
	
	// for each tile 
	for (int nb_tile = 0; nb_tile < tiles_half; nb_tile++) {

		// height of the tile
		for (int i = 0; i < TILE; i++) {
			height_index = (int)(nb_tile/width_tile)*TILE+i;
			// current imaginary value
			y = y_start + height_index * dy; 

			// width of the tile
			#pragma omp simd simdlen(64)
			for (int j = 0; j < TILE; j++) {
				width_index = (nb_tile%width_tile)*TILE+j;
				// current real value
				x = x_start + width_index * dx; 
				zReal = x;
				zImag = y;
				value = limit;
				for (int l = 0; l < limit; ++l) {
					r2 = zReal * zReal;
					i2 = zImag * zImag;
					zImag = 2.0f * zReal * zImag + y;
					zReal = r2 - i2 + x;
					if (r2 + i2 > 4.0f) {
						value = l;
						break;
					}
				}
				// update data with the computed value
				data[height_index*width+width_index] = value;
				data[(height-height_index-1)*width+width_index] = value;
			}
		}
	}
	return data;
}