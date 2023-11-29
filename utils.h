#pragma once

#include <stdlib.h>
#include <cstdio>

typedef unsigned int uint;

__attribute__((always_inline)) inline
static double linterp(double x, double xmin, double xmax, double a, double b) { return a + x/(xmax-xmin)*(b-a); }

__attribute__((always_inline)) inline
static int clampi(int x, int a, int b)            { if (x < a) { return a; } else if (x > b) { return b; } else { return x; } }

__attribute__((always_inline)) inline
static double clamp(double x, double a, double b) { if (x < a) { return a; } else if (x > b) { return b; } else { return x; } }

static void writearr(FILE* fp, double buf[], uint n) {
  for (uint i = 0; i < n; ++i) {
    // todo: a better way to dump the array than just the string
    fprintf(fp, "%.17f ", buf[i]);
  }
  fprintf(fp, "\n");
}

// static void linspace(double buf[], uint n, double xmin, double xmax) {
// 	// @brief fills a grid with a constant spacing from xmin to xmax
// 	double d = (xmax-xmin)/(n-1.0);
// 	for (size_t i = 0; i < n; ++i) {
// 		buf[i] = xmin + ((double) i)*d;
// 	}
// }

// static void full(double buf[], uint n, double val) {
// 	// @brief fills the buffer with n copies of val
// 	for (size_t i = 0; i < n; ++i) {
// 		buf[i] = val;
// 	}
// }

// static void map(double buf[], uint n, double (*fn)(double)) {
// 	for (size_t i = 0; i < n; ++i) {
// 		buf[i] = fn(buf[i]);
// 	}
// }

// static double reduce(double buf[], uint n, double (*fn)(double, double)) {
// 	double accum = 0.0;
// 	for (size_t i = 0; i < n; ++i) {
// 		accum = fn(accum, buf[i]);
// 	}
// 	return accum;
// }

