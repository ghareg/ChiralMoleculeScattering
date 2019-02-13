#include "kzIntegral.h"
#include "rIntegral.h"
#include <iostream>
#include <gsl/gsl_sf_hyperg.h>

int main()
{
	PauliMatrix pal;
	initPauliMatrix(pal);

	LValStruct LVal;
	GValStruct GVal;
	PGValStruct PGVal;
	Param pm;

	double E = 0.5;
	double k = sqrt(E);
	double alpha = 0.0;
	double beta = 0.0;
	double theta = 1.2;
	double tau = 1.2;
	Complex knm = 0.0;
	double qnm = 0.0;
	double absm = 0;
	double kSinthSq = 0.0;

	double fact[nmax + mmax + 1];
	for (int k = 0; k < nmax + mmax + 1; ++k) {
		if (k == 0) {
			fact[k] = 1.0;
		}
		else {
			fact[k] = k * fact[k - 1];
		}
	}

	pm.updateValues(k, alpha, beta, theta, tau);
	Matrix2cd ft;
	Matrix2cd f0;
	Matrix2cd f1;
	Matrix2cd f2;
	Matrix2cd f3;
	Matrix2cd f4;
	Matrix2cd rho = 0.5 * pal.pal0;
	Complex Pol;
	Complex I10t;
	Complex I21t;
	Complex I31t;
	Complex I23t;
	Complex I33t;
	Complex I42t;
	Complex I52t;
	Complex I2m1t;
	Complex I3m1t;
	Complex I40t;
	Complex I50t;
	Complex mult;
	for (int nc = 0; nc <= nmax; ++nc) {
		for (int mc = -mmax; mc <= mmax; ++mc) {
			absm = mc >=0 ? mc : -mc;
			qnm = 4.0 * nc + 2.0 * absm + 2.0;
			knm = std::sqrt(k * k - qnm);
			updateLVal(LVal, nc, mc, pm, pal);
			updateGVal(GVal, LVal, pm.kCosal, knm);
			updatePGVal(PGVal, GVal, pm.kCosth, pm.kCosal, knm);

			Complex I30i[nc + 1];
			Complex I31i[nc + 1];
			kSinthSq = pm.kSinth * pm.kSinth;
			I30i[0] = -pm.kSinth * gsl_sf_hyperg_1F1_int (1, (int) absm + 2, 0.5 * pm.kSinth * pm. kSinth) / (2.0 * (absm + 1));
			I31i[0] = -1.0 / pm.kSinth;
			for (int k = 1; k < nc + 1; ++k) {
				I30i[k] = -(powm1(k) / k) * pm.kSinth * Laguerre(k - 1, (int) absm + 1, kSinthSq) + ((k + absm) / k) * I30i[k - 1];
				I31i[k] = -(powm1(k) / pm.kSinth) * Laguerre(k, (int) absm - 1, kSinthSq) + I31i[k - 1];
			}
			if (std::abs(pm.kSinth) < thresh) {
				mult = std::pow(I, absm) * std::sqrt(fact[nc] * Pi / fact[nc + (int) absm]);
			}
			else {
				mult = std::pow(-pm.kSinth * I, absm) * std::sqrt(fact[nc] * Pi / fact[nc + (int) absm]) * std::exp(-0.5 * kSinthSq);
			}

			I10t = mult * I10(nc, -mc, -pm.kSinth, pm.tau);
			I21t = mult * I21(nc, -mc, -pm.kSinth, pm.tau);
			I31t = mult * I31(nc, -mc, -pm.kSinth, pm.tau);
			I33t = mult * I33(nc, -mc, -pm.kSinth, pm.tau);
			I42t = mult * I42(nc, -mc, -pm.kSinth, pm.tau);
			I52t = mult * I52(nc, -mc, -pm.kSinth, pm.tau);
			I2m1t = mult * I2m1(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
			I3m1t = mult * I3m1(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
			I40t = mult * I40(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
			I50t = mult * I50(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);

			f0 += I10t * Iz2(pm.kCosal, pm.kCosth, 0, Lz) * LVal.L1 + I10t * 
				Iz1(pm.kCosal, pm.kCosth, 0, Lz) * LVal.L2;
			f1 += F1 * I10t * PGVal.PG2 + (F2 * I21t + F3 * I31t) * PGVal.PG1;
			f2 += (2.0 * I31t + F3 * I10t) * PGVal.PG3 + (-I33t + absm * I31t + mc * 1.0 * I * I21t +
				I42t) * PGVal.PG4;
			f3 += (2.0 * I21t + F2 * I10t) * PGVal.PG3 + (-I23t + absm * I21t - mc * 1.0 * I * I31t +
				I52t) * PGVal.PG4;
			f4 += (mc * 2.0 * I * I10t - F2 * I31t + F2 * absm * I3m1t + mc * F2 * I * I2m1t + 
				F2 * I40t + F3 * I21t - F3 * absm * I2m1t + mc * F3 * I * I3m1t - F3 * I50t) * PGVal.PG1;
		}
	}
	f0 = (-1.0  / (4 * Pi)) * f0;
	f1 = (-1.0 / (4 * Pi)) * f1;
	f2 = (I * alpha / (4 * Pi)) * pal.pal1 * f2;
	f3 = (-I * alpha / (4 * Pi)) * pal.pal2 * f3;
	f4 = (I * alpha / (4 * Pi)) * pal.pal3 * f4;
	ft = f0 + f1 + f2 + f3 + f4;
	ft = ft * rho * ft.adjoint();
	Pol = (pal.pal3 * ft).trace() / ft.trace();
	std::cout << Pol.real() << "\t" << Pol.imag() << std::endl;
}
