#include "kzIntegral.h"
#include "rIntegral.h"

Complex Iz1(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		return I * (std::exp(I * (k1 - k2) * a) - std::exp(I * (k1 - k2) * b)) / (k1 - k2);
	}
	else {
		return (b - a);
	}
}

Complex Iz2(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		return ((1.0 - I * (k1 - k2) * b) * std::exp(I * (k1 - k2) * b) - (1.0 - I * (k1 - k2) * a) * std::exp(I * (k1 - k2) * a)) / 
			((k1 - k2) * (k1 - k2));
	}
	else {
		return 0.5 * (b * b - a * a);
	}
}

Complex Iz3(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		Complex mult1 = 2.0 * I + 2.0 * (k1 - k2) * b - I * (k1 -k2) * (k1 - k2) * b * b;
		Complex mult2 = 2.0 * I + 2.0 * (k1 - k2) * a - I * (k1 -k2) * (k1 - k2) * a * a;
		return (mult1 * std::exp(I * (k1 - k2) * b) - mult2 * std::exp(I * (k1 - k2) * a)) / 
			((k1 - k2) * (k1 - k2) * (k1 - k2));
	}
	else {
		return (b * b * b - a * a * a) / 3.0;
	}
}

void updateGVal(GValStruct& GVal, const LValStruct& LVal, double kCosal, Complex knm)
{
	if (std::abs(kCosal - knm) < thresh) {
	}
	else if (std::abs(kCosal + knm) < thresh) {
	}
	else {
		GVal.G1 = (I / (2.0 * knm * (kCosal - knm) * (kCosal - knm))) * LVal.L1 + 
			(1.0 / (2.0 * knm * (kCosal - knm))) * LVal.L2;
		GVal.G2 = (-2.0 * I * kCosal / ((kCosal * kCosal - knm * knm) * (kCosal * kCosal - knm * knm))) * LVal.L1 -
		   (1.0 / (kCosal * kCosal - knm * knm)) * LVal.L2;
		GVal.G3 = (-I * std::exp(I * (kCosal + knm) * Lz) / (2.0 * knm * (kCosal + knm) * (kCosal + knm)) - 
			Lz * std::exp(I * (kCosal + knm) * Lz) / (2.0 * knm * (kCosal + knm))) * LVal.L1 - 
			(std::exp(I * (kCosal + knm) * Lz) / (2.0 * knm * (kCosal + knm))) * LVal.L2;
		GVal.G4 = (-1.0 / (kCosal * kCosal - knm * knm)) * LVal.L1;
	}
}

void updatePGVal(PGValStruct& PGVal, const GValStruct& GVal, double kCosth, double kCosal, Complex knm)
{
	if (std::abs(kCosal - knm) < thresh) {
	}
	else if (std::abs(kCosal + knm) < thresh) {
	}
	else {
		Complex Iz1knmth = Iz1(knm, kCosth, 0, Lz);
		Complex Iz1alth = Iz1(kCosal, kCosth, 0, Lz);
		Complex Iz1mknmth = Iz1(-knm, kCosth, 0, Lz);
		Complex Iz2alth = Iz2(kCosal, kCosth, 0, Lz);
		Complex Iz2knmth = Iz2(knm, kCosth, 0, Lz);
		Complex Iz2mknmth = Iz2(-knm, kCosth, 0, Lz);
		Complex Iz3alth = Iz3(kCosal, kCosth, 0, Lz);
		PGVal.PG1 = GVal.G1 * Iz1knmth + GVal.G2 * Iz1alth +
			GVal.G3 * Iz1mknmth + GVal.G4 * Iz2alth;
		PGVal.PG2 = GVal.G1 * Iz2knmth + GVal.G2 * Iz2alth +
			GVal.G3 * Iz2mknmth + GVal.G4 * Iz3alth;
		PGVal.PG3 = I * knm * Iz1knmth * GVal.G1 + I * kCosal * Iz1alth * GVal.G2 -
			I * knm * Iz1mknmth * GVal.G3 + (Iz1alth + I * kCosal * Iz2alth) * GVal.G4;
		PGVal.PG4 = -I * kCosth * PGVal.PG1 + PGVal.PG3;
		PGVal.PG5 = -I * kCosth * PGVal.PG2 + I * knm * Iz2knmth * GVal.G1 + I * kCosal * Iz2alth * GVal.G2 +
			- I * knm * Iz2mknmth * GVal.G3 + (Iz2alth + I * kCosal * Iz3alth) * GVal.G4;
	}
}
