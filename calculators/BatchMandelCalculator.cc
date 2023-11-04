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
	tmp_data = (int *)(_mm_malloc(64 * sizeof(int), 64));
	// x_init = (float *)(_mm_malloc(64 * sizeof(float), 64));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	data = NULL;
	_mm_free(complex_tmp);
	complex_tmp = NULL;
	_mm_free(tmp_data);
	complex_tmp = NULL;
	// _mm_free(x_init);
	// x_init = NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {
	float *pcomplexReal;
	float *pcomplexImag;
	int *ptmp = tmp_data;
	// float *pxInit = x_init;
	float y, x, r2, i2;
	int sum;
	int TILE = 64;
	int TILE_sizeofint = TILE*sizeof(int);

	int N_width = width/TILE;

	// height loop
	for (int i = 0; i < height/2; i++) {

		// current imaginary value
		y = y_start + i * dy;

		
		// #pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
			for (int j = 0; j < N_width; j++) {
				
				// set the initial imaginary values
				std::uninitialized_fill(complex_tmp, complex_tmp + 2*TILE, y);
				// set the initial real values
				pcomplexReal = complex_tmp;
				// pxInit = x_init;
				int j_TILE = j*TILE;
				#pragma omp simd simdlen(64)
				for (int jj = 0; jj < TILE; jj++) {
					*pcomplexReal = x_start + (j_TILE+jj) * dx;
					// *(pxInit++) = *pcomplexReal;
					pcomplexReal+=2;
				}

				// set the initial tmp_data with limit
				std::uninitialized_fill(tmp_data, tmp_data + TILE, limit);


				// iterations
				for (int l = 0; l < limit; l++) {
					pcomplexReal = complex_tmp;
					pcomplexImag = complex_tmp+1;
					ptmp = tmp_data;
					// pxInit = x_init;
					sum = 0;

					// width loop 
					#pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
					// for (int j = 0; j < width/64; j++) {
						// inner loop with TILE size
						for (int jj = 0; jj < TILE; jj++) {
							// calculate one iteration of mandelbrot
							r2 = *pcomplexReal * *pcomplexReal;
							i2 = *pcomplexImag * *pcomplexImag;

							if (r2 + i2 > 4.0f) {
								if (*ptmp == limit)
									*ptmp = l;
								sum++;
							}

							// update complex tmp
							*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
							*pcomplexReal = r2 - i2 + (x_start + (j_TILE+jj) * dx);
				
							pcomplexImag+=2;
							pcomplexReal+=2;
							ptmp++;
							// pxInit++;
						}
						// end the iterations before reaching the limit
						if (sum == width) break;
					}
					memcpy(data + (i*width)+j_TILE, tmp_data, TILE_sizeofint);
					memcpy(data + (height-i-1)*width+j_TILE, tmp_data, TILE_sizeofint);
				}
	}

	return data;

	/*
	float *pcomplexReal = complex_tmp;
	float *pcomplexImag = complex_tmp+1;
	float y, x, r2, i2;
	int sum, height_offset, height_offset_reverse;
	int TILE = 64;
	int N_width = width/TILE;
	int N_height = height/(TILE*2);
	int height_index, width_index, value;

	std::uninitialized_fill(data, data + height*width, limit);

	for (int nb_tile = 0; nb_tile < ((width*height)/(TILE*TILE))/2; nb_tile++) {

		// initialize complex_tmp
		for (int i = 0; i < TILE; i++) {
			// #pragma omp simd simdlen(64)
			for (int j = 0; j < TILE; j++) {
				height_index = (int)(nb_tile/N_width)*TILE + i;
				width_index = (nb_tile%N_width)*TILE + j;

				complex_tmp[2*(i*TILE+j)] = y_start + height_index * dy;
				complex_tmp[2*(i*TILE+j)+1] = x_start + width_index * dx;
			}
		}

		sum = 0;

		// limit
		for (int l = 0; l < limit; l++) {

			// height
			for (int i = 0; i < TILE; i++) {
				height_index = (int)(nb_tile/N_width)*TILE + i;
				y = y_start + height_index * dy;

				// width
				// #pragma omp simd simdlen(64)
				for (int j = 0; j < TILE; j++) {
					width_index = (nb_tile%N_width)*TILE + j;

					// if (data[height_index*width+width_index] != limit) {
					// 	sum++;
					// 	continue;
					// }

					*pcomplexReal = complex_tmp[2*(i*TILE+j)];
					*pcomplexImag = complex_tmp[2*(i*TILE+j)+1];

					r2 = *pcomplexReal * *pcomplexReal;
					i2 = *pcomplexImag * *pcomplexImag;

					if (r2 + i2 > 4.0f) {
						if (data[height_index*width+width_index] == limit) {
							data[height_index*width+width_index] = l;
						}
						sum++;
					}

					complex_tmp[2*(i*TILE+j)+1] = 2.0f * *pcomplexReal * *pcomplexImag + y;
					x = x_start + width_index * dx;
					complex_tmp[2*(i*TILE+j)] = r2 - i2 + x;
				}	// width
			}	// height
		}
	}

	return data;

	*/


/*
		for (int ii = 0; ii < TILE; ii++) {
			// current imaginary value
			height_index = i * TILE + ii;
			height_offset = height_index * width;
			height_offset_reverse = (height - height_index - 1) * width;
			y = y_start + height_index * dy;

			

			for (int j = 0; j < N_width; j++) {

				// set the initial imaginary values
				std::uninitialized_fill(complex_tmp, complex_tmp + 2*TILE, y);
				// set the initial real values
				pcomplexReal = complex_tmp;
				// #pragma omp simd simdlen(64)
				for (int jj = 0; jj < TILE; jj++) {
					*pcomplexReal = x_start + (j * TILE + jj) * dx;
					pcomplexReal+=2;
				}

				// set the initial tmp_data with limit
				std::uninitialized_fill(tmp_data, tmp_data + TILE, limit);

				for (int l = 0; l < limit; l++) {
					pcomplexReal = complex_tmp;
					pcomplexImag = complex_tmp+1;
					ptmp = tmp_data;
					sum = 0;
					#pragma omp simd simdlen(64)
					for (int jj = 0; jj < TILE; jj++) {
						// calculate one iteration of mandelbrot
						r2 = *pcomplexReal * *pcomplexReal;
						i2 = *pcomplexImag * *pcomplexImag;

						if (r2 + i2 > 4.0f) {
							if (*ptmp == limit)
								*ptmp = l;
							sum++;
						}

						// update complex tmp
						*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
						width_index = j * TILE + jj;
						x = x_start + width_index * dx;
						*pcomplexReal = r2 - i2 + x;
			
						pcomplexImag+=2;
						pcomplexReal+=2;
						ptmp++;
					}
					if (sum == width) break;
				}
				// spravne nakopirovat data na spravne misto
				memcpy(data + height_offset + width_index, tmp_data, TILE*sizeof(int));
				memcpy(data + height_offset_reverse + width_index, tmp_data, TILE*sizeof(int));
			}
		}
	}

	return data;
	*/

	// 	// current imaginary value
	// 	y = y_start + i * dy;

	// 	// set the initial imaginary values
	// 	std::uninitialized_fill(complex_tmp, complex_tmp + 2*width, y);
	// 	// set the initial real values
	// 	pcomplexReal = complex_tmp;
	// 	#pragma omp simd simdlen(64)
	// 	for (int j = 0; j < width; j++) {
	// 		*pcomplexReal = x_start + j * dx;
	// 		pcomplexReal+=2;
	// 	}

	// 	// set the initial tmp_data with limit
	// 	std::uninitialized_fill(tmp_data, tmp_data + width, limit);
		
	// 	// iterations
	// 	for (int l = 0; l < limit; l++) {
	// 		pcomplexReal = complex_tmp;
	// 		pcomplexImag = complex_tmp+1;
	// 		ptmp = tmp_data;
	// 		sum = 0;

	// 		// width loop 
	// 		#pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
	// 		for (int j = 0; j < width; j++) {
	// 			// calculate one iteration of mandelbrot
	// 			r2 = *pcomplexReal * *pcomplexReal;
	// 			i2 = *pcomplexImag * *pcomplexImag;

	// 			if (r2 + i2 > 4.0f) {
	// 				if (*ptmp == limit)
	// 					*ptmp = l;
	// 				sum++;
	// 			}

	// 			// update complex tmp
	// 			*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
	// 			x = x_start + j * dx;
	// 			*pcomplexReal = r2 - i2 + x;
	
	// 			pcomplexImag+=2;
	// 			pcomplexReal+=2;
	// 			ptmp++;
	// 		}

	// 		// end the iterations before reaching the limit
	// 		if (sum == width) break;
	// 	}

	// 	memcpy(data + (i*width), tmp_data, width*sizeof(int));
	// 	memcpy(data + (height-i-1)*width, tmp_data, width*sizeof(int));
	// }

	// return data;


	// // height loop
	// for (int i = 0; i < N_height; i++) {
	// 	for (int ii = 0; ii < TILE; ii++) {
	// 		// current imaginary value
	// 		height_index = i * TILE + ii;
	// 		height_offset = height_index * width;
	// 		height_offset_reverse = (height - height_index - 1) * width;
	// 		y = y_start + height_index * dy;

	// 		// set the initial imaginary values
	// 		std::uninitialized_fill(complex_tmp, complex_tmp + 2*TILE, y);
	// 		// set the initial real values
	// 		pcomplexReal = complex_tmp;
	// 		#pragma omp simd simdlen(64)
	// 		for (int j = 0; j < N_width; j++) {
	// 			*pcomplexReal = x_start + j * dx;
	// 			pcomplexReal+=2;
	// 		}

	// 		// set the initial tmp_data with limit
	// 		std::uninitialized_fill(tmp_data, tmp_data + TILE, limit);
			
	// 		// iterations
	// 		for (int l = 0; l < limit; l++) {
	// 			pcomplexReal = complex_tmp;
	// 			pcomplexImag = complex_tmp+1;
	// 			ptmp = tmp_data;
	// 			sum = 0;

	// 			// width loop 
	// 			for (int j = 0; j < N_width; j++) {
	// 				ptmp = tmp_data;
	// 				#pragma omp simd simdlen(64) aligned(pcomplexReal, pcomplexImag: 64) reduction(+:sum)
	// 				for (int jj = 0; jj < TILE; jj++) {
	// 					width_index = j * TILE + jj;
	// 					x = x_start + width_index * dx;
				
	// 					// calculate one iteration of mandelbrot
	// 					r2 = *pcomplexReal * *pcomplexReal;
	// 					i2 = *pcomplexImag * *pcomplexImag;

	// 					if (r2 + i2 > 4.0f) {
	// 						if (*ptmp == limit)
	// 							*ptmp = l;
	// 						sum++;
	// 					}

	// 					// update complex tmp
	// 					*pcomplexImag = 2.0f * *pcomplexReal * *pcomplexImag + y;
	// 					*pcomplexReal = r2 - i2 + x;
	
	// 					pcomplexImag+=2;
	// 					pcomplexReal+=2;
	// 					ptmp++;
	// 				}
	// 			}
				

	// 			// end the iterations before reaching the limit
	// 			if (sum == width) break;
	// 		}
	// 			memcpy(data + (height_offset), tmp_data, TILE*sizeof(int));
	// 			memcpy(data + (height_offset_reverse)*width, tmp_data, TILE*sizeof(int));
			
	// 	}
	// }

	// return data;




/*
	
	int *pdata = data;

	int TILE = 64;	// TODO
	int N_width = width/TILE;
	int N_height = height/(TILE*2);
	int height_index, width_index, value;
	float x, y, zReal, zImag, r2, i2;
	int height_offset, height_offset_reverse;

	// height
	for (int i = 0; i < N_height; i++) {
		for (int ii = 0; ii < TILE; ii++) {
			height_index = i * TILE + ii;
			height_offset = height_index * width;
			height_offset_reverse = (height - height_index - 1) * width;
			y = y_start + height_index * dy;

			// width
			for (int j = 0; j < N_width; j++) {

				#pragma omp simd
				for (int jj = 0; jj < TILE; jj++) {
					width_index = j * TILE + jj;
					x = x_start + width_index * dx;

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
					pdata[height_offset + width_index] = value;
					pdata[height_offset_reverse + width_index] = value;
				}
			}
		}
	}
	return data;
	*/
}


