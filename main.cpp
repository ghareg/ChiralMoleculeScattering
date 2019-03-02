#include "kzIntegral.h"
#include "rIntegral.h"
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_hyperg.h>
#include <thread>
void calcFullPol(double* Pol, const double* fact, const PauliMatrix& pal);
void calcPol(double* Pol, const double* fact, const PauliMatrix& pal, int rs, int re);

int main()
{
	PauliMatrix pal;
	initPauliMatrix(pal);

	double* Pol = new double[ND * ND]; 
	double fact[nmax + mmax + 1];
	for (int k = 0; k < nmax + mmax + 1; ++k) {
		if (k == 0) {
			fact[k] = 1.0;
		}
		else {
			fact[k] = k * fact[k - 1];
		}
	}

	calcFullPol(Pol, fact, pal);
	double E = 0.001;
	double theta = 0.0;
	std::cout << std::setprecision(4);
	for (int ind1 = 0; ind1 < ND; ++ind1) {
		E = 0.001;
		for (int ind2 = 0; ind2 < ND; ++ind2) {
			std::cout << theta << "\t" << E << "\t" << Pol[ind2 * ND + ind1] << std::endl;
			E += EMax / ND;
		}
		std::cout << std::endl;
		theta += Pi / ND;
	}
	delete[] Pol;
}

void calcFullPol(double* Pol, const double* fact, const PauliMatrix& pal)
{
	double* PolMat[nCore];
	int div = ND / nCore;
	int size1 = div * ND;
	int size2 = (ND - (nCore - 1) * div) * ND;
	for (int i = 0; i < nCore - 1; ++i) {
		PolMat[i] = new double[size1];
	}
	PolMat[nCore - 1] = new double[size2];
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcPol, PolMat[i], fact, std::cref(pal), i * div, (i + 1) * div);
	}
						
	calcPol(PolMat[nCore - 1], fact, pal, (nCore - 1) * div, ND);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}

	int ind = 0;
	for (int i = 0; i < nCore - 1; ++i) {
		for (int j = 0; j < size1; ++j) {
			Pol[ind++] = PolMat[i][j];
		}
	}
	for (int j = 0; j < size2; ++j) {
		Pol[ind++] = PolMat[nCore - 1][j];
	}
	for (int i = 0; i < nCore; ++i) {
		delete[] PolMat[i];
	}
}

void calcPol(double* Pol, const double* fact, const PauliMatrix& pal, int rs, int re)
{
	double E = EnMult * (rs * EMax / ND + 0.001);
	double k = sqrt(E);
	double alpha = 0.0;
	double beta = 0.0;
	double theta = 0;
	double tau = 0;
	Complex knm = 0.0;
	double qnm = 0.0;
	double absm = 0;
	double kSinthSq = 0.0;
	Matrix2cd ft;
	Matrix2cd f0;
	Matrix2cd f1;
	Matrix2cd f2;
	Matrix2cd f3;
	Matrix2cd f4;
	Matrix2cd rho = 0.5 * pal.pal0;
	Complex I10t;
	Complex I12t;
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
	Complex I60t;
	Complex I62t;
	Complex I70t;
	Complex I72t;
	Complex I451t;
	Complex I91t;
	Complex I101t;
	Complex mult;
	LValStruct LVal;
	GValStruct GVal;
	PGValStruct PGVal;
	Param pm;
	double Polc = 0.0;

	for (int ind1 = 0; ind1 < (re - rs); ++ind1) {
		k = sqrt(E);
		theta = 0.0;
		for (int ind2 = 0; ind2 < ND; ++ind2) {
			tau = 0.0;
			Polc = 0.0;	
			for (int ind3 = 0; ind3 < ND; ++ind3) {
				beta = 0.0;
				for (int ind4 = 0; ind4 < ND; ++ind4) {
					pm.updateValues(k, alpha, beta, theta, tau);
					f0 << 0.0, 0.0, 0.0, 0.0;
					f1 << 0.0, 0.0, 0.0, 0.0;
					f2 << 0.0, 0.0, 0.0, 0.0;
					f3 << 0.0, 0.0, 0.0, 0.0;
					f4 << 0.0, 0.0, 0.0, 0.0;
					for (int nc = 0; nc <= nmax; ++nc) {
						for (int mc = -mmax; mc <= mmax; ++mc) {
							absm = mc >=0 ? mc : -mc;
							qnm = 4.0 * nc + 2.0 * absm + 2.0;
							knm = std::sqrt(Complex(k * k - qnm, 0.0));
							updateLVal(LVal, nc, mc, pm, pal, fact);
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
							I12t = mult * I12(nc, -mc, -pm.kSinth, pm.tau);
							I21t = mult * I21(nc, -mc, -pm.kSinth, pm.tau);
							I31t = mult * I31(nc, -mc, -pm.kSinth, pm.tau);
							I23t = mult * I23(nc, -mc, -pm.kSinth, pm.tau);
							I33t = mult * I33(nc, -mc, -pm.kSinth, pm.tau);
							I42t = mult * I42(nc, -mc, -pm.kSinth, pm.tau);
							I52t = mult * I52(nc, -mc, -pm.kSinth, pm.tau);
							I2m1t = mult * I2m1(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I3m1t = mult * I3m1(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I40t = mult * I40(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I50t = mult * I50(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I60t = mult * I60(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I62t = mult * I62(nc, -mc, -pm.kSinth, pm.tau);
							I70t = mult * I70(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I72t = mult * I72(nc, -mc, -pm.kSinth, pm.tau);
							I451t = mult * I451(nc, -mc, -pm.kSinth, pm.tau);
							I91t = mult * I91(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);
							I101t = mult * I101(nc, -mc, -pm.kSinth, pm.tau, I30i, I31i);	


							f0 += I10t * Iz2(pm.kCosal, pm.kCosth, 0, Lz) * LVal.L1 + I10t * 
								Iz1(pm.kCosal, pm.kCosth, 0, Lz) * LVal.L2;
							f1 += (1.0 + I * kappaD) * F1 * I10t * PGVal.PG2 + (-I * kappa * I12t + (1.0 + I * kappaD) * (F2 * I21t + F3 * I31t)) * PGVal.PG1;
							f2 += (2.0 * (1.0 + I * kappa) * I31t + F3 * I10t) * PGVal.PG3 + ((1.0 + I * kappa) * (-I33t + absm * I31t + mc * 1.0 * I * I21t +
								I42t) + F2 * (-0.5 * I62t + 0.5 * absm * I60t + 0.5 * mc * I * (I10t + I70t) + 0.5 * I91t) +
								F3 * (-0.5 * (I12t - I72t) + 0.5 * absm * (I10t - I70t) + 0.5 * mc * I * I60t + 0.5 * (I451t - I101t))) * PGVal.PG4
								+ F1 * (-I31t + absm * I3m1t + 1.0 * mc * I * I2m1t + I40t) * PGVal.PG5;
							f3 += (2.0 * (1.0 + I * kappa) * I21t + F2 * I10t) * PGVal.PG3 + ((1.0 + I * kappa) * (-I23t + absm * I21t - mc * 1.0 * I * I31t +
								I52t) + F2 * (-0.5 * (I12t + I72t) + 0.5 * absm * (I10t + I70t) - 0.5 * mc * I * I60t + 0.5 * (I101t + I451t)) +
								F3 * (-0.5 * I62t + 0.5 * absm * I60t - 0.5 * mc * I * (I10t - I70t) + 0.5 * I91t)) * PGVal.PG4
								+ F1 * (-I21t + absm * I2m1t - 1.0 * mc * I * I3m1t + I50t) * PGVal.PG5;
							f4 += (mc * 2.0 * (1.0 + I * kappa) * I * I10t - F2 * I31t + F2 * absm * I3m1t + mc * F2 * I * I2m1t + 
								F2 * I40t + F3 * I21t - F3 * absm * I2m1t + mc * F3 * I * I3m1t - F3 * I50t) * PGVal.PG1;
						}
					}
					f0 = (-1.0  / (4 * Pi)) * f0;
					f1 = (-1.0 / (4 * Pi)) * f1;
					f2 = (alphaSOC / (4 * Pi)) * I * (pal.pal1 * f2);
					f3 = (-alphaSOC / (4 * Pi)) * I * pal.pal2 * f3;
					f4 = (alphaSOC / (4 * Pi)) * I * pal.pal3 * f4;
					ft = f0 + f1 + f2 + f3 + f4;
					ft = ft * rho * ft.adjoint();
					Polc += ((pal.pal3 * ft).trace() / ft.trace()).real();
					beta += 2 * Pi / ND;
				}
				tau += 2 * Pi / ND;
			}
			Pol[ind1 * ND + ind2] = Polc / (1.0 * ND * ND);
			theta += Pi / ND;
		}
		E += EnMult * EMax / ND;
	}
}
