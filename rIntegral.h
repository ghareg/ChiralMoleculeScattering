#ifndef RINTEGRAL_H_
#define RINTEGRAL_H_

#include "LDefinit.h"

double powm1(int n);
double Laguerre(int n, int m, double val);

Complex I10(int n, int m, double kSinal, double beta);
Complex I12(int n, int m, double kSinal, double beta);
Complex I21(int n, int m, double kSinal, double beta);
Complex I23(int n, int m, double kSinal, double beta);
Complex I2m1(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31);
Complex I31(int n, int m, double kSinal, double beta);
Complex I33(int n, int m, double kSinal, double beta);
Complex I3m1(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31);
Complex I40(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31);
Complex I42(int n, int m, double kSinal, double beta);
Complex I50(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31);
Complex I52(int n, int m, double kSinal, double beta);
Complex I60(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31);
#endif
