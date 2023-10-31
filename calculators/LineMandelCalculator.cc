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

	for (int i = 0; i < height/2; ++i) {
		y = y_start + i * dy; // current imaginary value
		// pdata = &data[i*width];
		
		for (int l = 0; l < limit; ++l) {
			pdata = &data[i*width];
			// compute iterations for one row
			// bool b = true;
			// #pragma omp simd
			// #pragma omp simd reduction (|:b)
			for (int j = 0; j < width; ++j) {
				// c uz ma dopocitanou max iteraci
				c[j] = (l == 0) ? false : c[j];
				
				// TODO b
				if (c[j]) {
					pdata++;
					continue;
				}

				// one pixel
				x = x_start + j * dx; // current real value

				zRealA[j] = (l == 0) ? x : zRealA[j];
				zImagA[j] = (l == 0) ? y : zImagA[j];

				r2 = zRealA[j] * zRealA[j];
				i2 = zImagA[j] * zImagA[j];

				*pdata = l;
				data[((height-i-1)*width)+j] = l;
				pdata++;

				zImagA[j] = 2.0f * zRealA[j] * zImagA[j] + y;
				zRealA[j] = r2 - i2 + x;
				
				if (r2 + i2 > 4.0f) {
					// b |= false;	
					c[j] = true;
				}
			}
			// if (!b) break;
		}
	}
	return data;
}


// int * LineMandelCalculator::calculateMandelbrot () {
// 	// @TODO implement the calculator & return array of integers
// 	int *pdata = data;
// 	float *pzReal = zReal;
// 	float *pzImag = zImag;

// 	// TODO remove
// 	for (int i = 0; i < height; i++) {
// 		#pragma omp simd
// 		for (int j = 0; j < width; j++){

// 			float x = x_start + j * dx; // current real value
// 			float y = y_start + i * dy; // current imaginary value
			
// 			*(pzReal) = x; // jen poprve
// 			*(pzImag) = y; // jen poprve
// 			pzReal++;
// 			pzImag++;
// 			*(pdata) = 0;
// 			pdata++;

// 		}
// 	}

// 	pdata = data;

// 	// return data;


// 	for (int iteration = 0; iteration < limit; ++iteration) {

// 		// std::cout << iteration << std::endl;

// 		for (int i = 0; i < height; i++) {

// 			#pragma omp simd 
// 			for (int j = 0; j < width; j++){

// 				float x = x_start + j * dx; // current real value
// 				float y = y_start + i * dy; // current imaginary value
				
// 				// mandelbrot - jedna iterace
// 				auto r2 = *pzReal * *pzReal;		// potreba uz nova hodnota
// 				auto i2 = *pzImag * *pzImag;		// potreba uz nova hodnota

// 				if (r2 + i2 < 4.0f)
// 					// return i;
// 					*(pdata) += 1;		// TODO add reduction na iterations (pozor, je to pole a ne jen promenna)
// 				// std::cout << "pdata: " << *pdata << std::endl;
// 				pdata++; 
				
// 				// ulozit do pole mezivysledku pro dalsi iteraci
// 				*pzImag = 2.0f * *pzReal * *pzImag + y;
// 				*pzReal = r2 - i2 + x;

// 				// std::cout << "iteration: " << iteration << ", x: " << x << ", y: " << y << std::endl;
// 				// std::cout << "pzReal: " << *pzReal << ", pzImag: " << *pzImag << std::endl;

// 				pzReal++;
// 				pzImag++;
// 			}
// 		}
// 		pdata = data;
// 		pzReal = zReal;
// 		pzImag = zImag;

// // for width 
// // for height 
// // spocti mandelbrota jedno opakovani
// // if (neco%na%druhou mensi 4.0)
// // 			sum+=1;
// // 			*(iterations++) += 1;		// TODO add reduction na iterations (pozor, je to pole a ne jen promenna)
		
// 	}

// 	return data;
// }


// int * LineMandelCalculator::calculateMandelbrot () {
// 	// @TODO implement the calculator & return array of integers
// 	int *pdata = data;

// 		#pragma omp simd reduction(+: cnt)
// 		for (int i = 0; i < height; i++)
// 		{
// 			#pragma omp simd
// 			for (int j = 0; j < width; j++)
// 			{
// 				float x = x_start + j * dx; // current real value
// 				float y = y_start + i * dy; // current imaginary value

// 				int value = mandelbrot(x, y, limit);

// 				// mandelbrot
// 				//T zReal = real;
// 				//T zImag = imag;
// 				auto r2 = x * x;		// potreba uz nova hodnota
// 				auto i2 = y * y;		// potreba uz nova hodnota

// 				if (r2 + i2 > 4.0f)
// 					return i;

				
// 				x = 2.0f * x * y + y;
// 				// zImag *= 2.0f * zReal;
// 				// zImag += imag;

// 				x = r2 - i2 + x;

// 				value = iteration;

// 				*(pdata++) = value;
// 				// pdata[]
// 			}
// 		}
// 	}
// 	return data;
// }

