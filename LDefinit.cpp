#include "rIntegral.h"
#include "LDefinit.h"

void initPauliMatrix(PauliMatrix& pal)
{
	pal.pal0 << 1.0, 0.0, 0.0, 1.0;
	pal.pal1 << 0.0, 1.0, 1.0, 0.0;
	pal.pal2 << 0.0, -I, I, 0.0;
	pal.pal3 << 1.0, 0.0, 0.0, -1.0;
}

void updateLVal(LValStruct& LVal, int n, int m, const Param& pm, const PauliMatrix& pal, const double* fact)
{
	Complex mult;
	int absm = m >=0 ? m : -m;
	if (std::abs(pm.kSinal) < thresh) {
		mult = std::pow(I, absm) * std::sqrt(fact[n] * Pi / fact[n + absm]);
	}
	else {
		mult = std::pow(-pm.kSinal * I, absm) * std::sqrt(fact[n] * Pi / fact[n + absm]) * std::exp(-0.5 * pm.kSinal * pm.kSinal);
	}
	Complex I1 = mult * I10(n, m, pm.kSinal, pm.beta); 
	LVal.L1 = F1 * I1 * pal.pal0;
	Complex I2 = mult * I21(n, m, pm.kSinal, pm.beta);
	Complex I3 = mult * I31(n, m, pm.kSinal, pm.beta);
	Complex I12 = mult * I12(n, m, pm.kSInal, pm.beta);
	LVal.L2 = (-I * kappa * I12 + F2 * I2 + F3 * I3) * pal.pal0 + alphaSOC * (2.0 * (1.0 + I * kappa) * pm.kCosal * I3 + (F3 * pm.kCosal - F1 * pm.kSinal * pm.sinBeta) * I1)  * pal.pal1 -
		alphaSOC * (2.0 * (1.0 + I * kappa) * pm.kCosal * I2 + (F2 * pm.kCosal - F1 * pm.kSinal * pm.cosBeta) * I1) * pal.pal2 + 
		alphaSOC * (2.0 * (1.0 + I * kappa) * pm.kSinal * (I2 * pm.sinBeta - I3 * pm.cosBeta) + (F2 * pm.kSinal * pm.sinBeta - F3 * pm.kSinal * pm.cosBeta) * I1) * pal.pal3;
}	

void Param::updateValues(double kc, double alphac, double betac, double thetac, double tauc)
{
	k = kc;
	alpha = alphac;
	beta = betac;
	theta = thetac;
	tau = tauc;
	kSinal = k * std::sin(alpha);
	kCosal = k * std::cos(alpha);
	kSinth = k * std::sin(theta);
	kCosth = k * std::cos(theta);
	sinBeta = std::sin(beta);
	cosBeta = std::cos(beta);
	sinTau = std::sin(tau);
	cosTau = std::cos(tau);
}
