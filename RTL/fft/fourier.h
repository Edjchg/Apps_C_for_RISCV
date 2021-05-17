  
#ifndef __FOURIER_H__
#define __FOURIER_H__


#include "complex.h"
#include <iostream>

void radix2DitCooleyTykeyFft(int k, int* indices, Complex* x, Complex* f) ;

#endif