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
	Complex I1 = mult * I10(n, m, pm.kSinal);
	Complex I1a1 = mult * I10a1(n, m, pm.kSinal);
	Complex I1s1 = mult * I10s1(n, m, pm.kSinal);
	LVal.L1 = F1 * I1 * pal.pal0;
	Complex I2 = mult * I21(n, m, pm.kSinal);
	Complex I2a1 = mult * I21a1(n, m, pm.kSinal);
	Complex I2s1 = mult * I21s1(n, m, pm.kSinal);
	Complex I3 = mult * I31(n, m, pm.kSinal);
	Complex I3a1 = mult * I31a1(n, m, pm.kSinal);
	Complex I3s1 = mult * I31s1(n, m, pm.kSinal);
	LVal.L2 = (F2 * I2 + F3 * I3) * pal.pal0 + alphaSOC * (2.0 * pm.kCosal * I3 + (F3 * pm.kCosal) * I1 + 0.5 * F1 * pm.kSinal * I * (I1a1 - I1s1)) * pal.pal1 -
		alphaSOC * (2.0 * pm.kCosal * I2 + (F2 * pm.kCosal) * I1 - 0.5 * F1 * pm.kSinal * (I1a1 + I1s1)) * pal.pal2 + 
		alphaSOC * (pm.kSinal * (-I * (I2a1 - I2s1) - I3a1 - I3s1) + (-0.5 * F2 * pm.kSinal * I * (I1a1 - I1s1) - 0.5 * F3 * pm.kSinal * (I1a1 + I1s1))) * pal.pal3;
}	

void Param::updateValues(double kc, double alphac, double thetac)
{
	k = kc;
	alpha = alphac;
	theta = thetac;
	kSinal = k * std::sin(alpha);
	kCosal = k * std::cos(alpha);
	kSinth = k * std::sin(theta);
	kCosth = k * std::cos(theta);
}
