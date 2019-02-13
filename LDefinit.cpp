#include "rIntegral.h"
#include "LDefinit.h"

void initPauliMatrix(PauliMatrix& pal)
{
	pal.pal0 << 1.0, 0.0, 0.0, 1.0;
	pal.pal1 << 0.0, 1.0, 1.0, 0.0;
	pal.pal2 << 0.0, -I, I, 0.0;
	pal.pal3 << 1.0, 0.0, 0.0, -1.0;
}

void updateLVal(LValStruct& LVal, int n, int m, const Param& pm, const PauliMatrix& pal)
{
	LVal.L1 = F1 * I10(n, m, pm.kSinal, pm.beta) * pal.pal0;
	Complex I2 = I21(n, m, pm.kSinal, pm.beta);
	Complex I3 = I31(n, m, pm.kSinal, pm.beta);	
	LVal.L2 = (F2 * I2 + F3 * I3) * pal.pal0 + 2.0 * alpha * pm.kCosal * I3 * pal.pal1 -
		2.0 * alpha * pm.kCosal * I2 * pal.pal2 + 2.0 * alpha * pm.kSinal * (I2 * pm.sinBeta - I3 * pm.cosBeta) * pal.pal3;
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
