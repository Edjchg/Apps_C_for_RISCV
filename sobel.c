/* Copyright 2019 SiFive, Inc */
/* SPDX-License-Identifier: Apache-2.0 */


#include <stdio.h>
#include <time.h>
#include "image.h"

#define	ABS(A)	((A)<(0) ? (-(A)):(A)) //ABS macro, faster than procedure

//Exact
void sobel(int* out, int a, int b, int c, int d, int f, int g, int h, int i) {

	int x1 = a + g; //8 -> 9
	int x2 = d * 2; //8 -> 9
	int x3 = x1 + x2; //9 -> 10

	int x4 = c + i; //8 -> 9
	int x5 = f * 2; //8 -> 9
	int x6 = x4 + x5; //9 -> 10

	int gx = x6 - x3; 
	int absGx = ABS(gx);

	int y1 = a + c; //8 -> 9
	int y2 = b * 2; //8 -> 9
	int y3 = y1 + y2; //9 -> 10

	int y4 = g + i; //8 -> 9
	int y5 = h * 2; //8 -> 9
	int y6 = y4 + y5; //9 -> 10

	int gy = y3 - y6;
	int absGy = ABS(gy);

	int res = absGx + absGy;

  if (res > 255) res = 255;
  
  *out = res;
}

//Approximate1, d = a
void sobelSW1(int* out, int a, int b, int c, int f, int g, int h, int i) 
{

	int x1 = a + g; //8 -> 9
	int x2 = a * 2; //8 -> 9 x2 = d * 2
	int x3 = x1 + x2; //9 -> 10

	int x4 = c + i; //8 -> 9
	int x5 = f * 2; //8 -> 9
	int x6 = x4 + x5; //9 -> 10

	int gx = x6 - x3; 
	int absGx = ABS(gx);

	int y1 = a + c; //8 -> 9
	int y2 = b * 2; //8 -> 9
	int y3 = y1 + y2; //9 -> 10

	int y4 = g + i; //8 -> 9
	int y5 = h * 2; //8 -> 9
	int y6 = y4 + y5; //9 -> 10

	int gy = y3 - y6;
	int absGy = ABS(gy);

	int res = absGx + absGy;

  if (res > 255) res = 255;
  
  *out = res;
}

//Approximate2, d = a, b = c
void sobelSW2(int* out, int a, int c, int f, int g, int h, int i) 
{

	int x1 = a + g; //8 -> 9
	int x2 = a * 2; //8 -> 9 x2 = d * 2
	int x3 = x1 + x2; //9 -> 10

	int x4 = c + i; //8 -> 9
	int x5 = f * 2; //8 -> 9
	int x6 = x4 + x5; //9 -> 10

	int gx = x6 - x3; 
	int absGx = ABS(gx);

	int y1 = a + c; //8 -> 9
	int y2 = c * 2; //8 -> 9 y2 = b * 2
	int y3 = y1 + y2; //9 -> 10

	int y4 = g + i; //8 -> 9
	int y5 = h * 2; //8 -> 9
	int y6 = y4 + y5; //9 -> 10

	int gy = y3 - y6;
	int absGy = ABS(gy);

	int res = absGx + absGy;

  if (res > 255) res = 255;
  
  *out = res;
}


//Approximate3, d = a, b = c, f = i
void sobelSW3(int* out, int a, int c, int g, int h, int i) 
{

	int x1 = a + g; //8 -> 9
	int x2 = a * 2; //8 -> 9 x2 = d * 2
	int x3 = x1 + x2; //9 -> 10

	int x4 = c + i; //8 -> 9
	int x5 = i * 2; //8 -> 9 x5 = f * 2
	int x6 = x4 + x5; //9 -> 10

	int gx = x6 - x3; 
	int absGx = ABS(gx);

	int y1 = a + c; //8 -> 9
	int y2 = c * 2; //8 -> 9 y2 = b * 2
	int y3 = y1 + y2; //9 -> 10

	int y4 = g + i; //8 -> 9
	int y5 = h * 2; //8 -> 9
	int y6 = y4 + y5; //9 -> 10

	int gy = y3 - y6;
	int absGy = ABS(gy);

	int res = absGx + absGy;

  if (res > 255) res = 255;
  
  *out = res;
}


//Approximate4, d = a, b = c, f = i, h = g
void sobelSW4(int* out, int a, int c, int g, int i) 
{

	int x1 = a + g; //8 -> 9
	int x2 = a * 2; //8 -> 9 x2 = d * 2
	int x3 = x1 + x2; //9 -> 10

	int x4 = c + i; //8 -> 9
	int x5 = i * 2; //8 -> 9 x5 = f * 2
	int x6 = x4 + x5; //9 -> 10

	int gx = x6 - x3; 
	int absGx = ABS(gx);

	int y1 = a + c; //8 -> 9
	int y2 = c * 2; //8 -> 9 y2 = b * 2
	int y3 = y1 + y2; //9 -> 10

	int y4 = g + i; //8 -> 9
	int y5 = g * 2; //8 -> 9 y5 = h * 2
	int y6 = y4 + y5; //9 -> 10

	int gy = y3 - y6;
	int absGy = ABS(gy);

	int res = absGx + absGy;

  if (res > 255) res = 255;
  
  *out = res;
}



int main() {

  int A = 0;
  int B = 0;
  int C = 0;
  int D = 0;
  int F = 0;
  int G = 0;
  int H = 0;
  int I = 0;
  int i,j;
  int out = 0;

  int RESULT[HEIGHT][WIDTH] = {0};

  //clock_t begin = clock();

  for (i = 1; i < HEIGHT-1; i++) {
    for (j = 1; j < WIDTH-1; j++) {
      A = IMAGE[i - 1][j - 1];
      B = IMAGE[i - 1][j];
      C = IMAGE[i - 1][j + 1];
      D = IMAGE[i][j-1];
      F = IMAGE[i][j + 1];
      G = IMAGE[i + 1][j - 1];
      H = IMAGE[i + 1][j];
      I = IMAGE[i + 1][j + 1];

      sobel(&out, A, B, C, D, F, G, H, I);

      RESULT[i][j] = out;
      printf("New pixel %i \n", out);
    }
  }

  //clock_t end = clock();

  //int clock_cycles = (int)(end - begin);
  //printf("clock cycles: %d\n",clock_cycles);

  printf("Sobel filtering\n");

	return 0;
}

