#ifndef RINTEGRAL_H_
#define RINTEGRAL_H_

#include "LDefinit.h"

double powm1(int n);
double Laguerre(int n, int m, double val);

Complex I10(int n, int m, double kSinal);
Complex I10a1(int n, int m, double kSinal);
Complex I10s1(int n, int m, double kSinal);
Complex I12(int n, int m, double kSinal);
Complex I21(int n, int m, double kSinal);
Complex I21a1(int n, int m, double kSinal);
Complex I21s1(int n, int m, double kSinal);
Complex I23(int n, int m, double kSinal);
Complex I2m1(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I31(int n, int m, double kSinal);
Complex I31a1(int n, int m, double kSinal);
Complex I31s1(int n, int m, double kSinal);
Complex I33(int n, int m, double kSinal);
Complex I3m1(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I40(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I42(int n, int m, double kSinal);
Complex I50(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I52(int n, int m, double kSinal);
Complex I60(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I62(int n, int m, double kSinal);
Complex I70(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I72(int n, int m, double kSinal);
Complex I451(int n, int m, double kSinal);
Complex I91(int n, int m, double kSinal, Complex* I30, Complex* I31);
Complex I101(int n, int m, double kSinal, Complex* I30, Complex* I31);
#endif
