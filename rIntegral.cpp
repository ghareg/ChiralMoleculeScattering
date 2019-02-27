#include "rIntegral.h"
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

double powm1(int n)
{
	if (n % 2) {
		return -1.0;
	}
	return 1.0;
}

double Laguerre(int n, int m, double val)
{
	if (n < 0) {
		return 0.0;
	}

	if (m > -1) {
		return  gsl_sf_laguerre_n(n, m, val);
	}
	else {
		return Laguerre(n, m + 1, val) - Laguerre(n - 1, m + 1, val);
	}
}

Complex I10(int n, int m, double kSinal, double beta)
{
	int absm = m >= 0 ? m : -m;
	if (std::abs(kSinal) < thresh) {
		if (m == 0) {
			return 2.0 * powm1(n);
		}
 		return 0.0;
	}
	
	return 2.0 * powm1(n) * std::exp(-I * (1.0 * m) * beta) * Laguerre(n, absm, kSinal * kSinal);
}

Complex I12(int n, int m, double kSinal, double beta)
{
	int absm = m >= 0 ? m : -m;
	if (std::abs(kSinal) < thresh) {
		if (m == 0) {
			return 4.0 * powm1(n) * (2 * n + 1);
		}
 		return 0.0;
	}
	
	return 2.0 * powm1(n) * std::exp(-I * (1.0 * m) * beta) * (-kSinal * kSinal * (Laguerre(n, absm + 2, kSinal * kSinal) +
				2.0 * Laguerre(n - 1, absm + 2, kSinal * kSinal) + Laguerre(n - 2, absm + 2, kSinal * kSinal)) +
			2.0 * (absm + 1) * (Laguerre(n, absm + 1, kSinal * kSinal) + Laguerre(n - 1, absm + 1, kSinal * kSinal)));
}

Complex I21(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			return -2.0 * (n + 1) * powm1(n) * I;
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	double kSinalSq = kSinal * kSinal;
	return powm1(n)* (kSinal * std::exp(- I * absm1 * beta) * (Laguerre(n, absm + 1, kSinalSq) +
				Laguerre(n - 1, absm + 1, kSinalSq)) - std::exp(- I * absm2 * beta) * 
			((n + absm) * Laguerre(n, absm - 1, kSinalSq) + (n + 1) * Laguerre(n + 1, absm - 1, kSinalSq))/kSinal) * I;
}

Complex I23(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			return 8.0 * (n + 1) * (n + 1) * powm1(n + 1) * I;
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	double kSinalSq = kSinal * kSinal;
	return powm1(n + 1) * I * (std::exp(- I * absm1 * beta) * (kSinalSq * kSinal * (Laguerre(n, absm + 3, kSinalSq) +
			3.0 * Laguerre(n - 1, absm + 3, kSinalSq) + 3.0 * Laguerre(n - 2, absm + 3, kSinalSq) + Laguerre(n - 3, absm + 3, kSinalSq)) 
			- 2.0 * kSinal * (absm + 2) * (Laguerre(n, absm + 2, kSinalSq) + 2.0 * Laguerre(n - 1, absm + 2, kSinalSq) + 
			Laguerre(n - 2, absm + 2, kSinalSq))) - std::exp(- I * absm2 * beta) * 
		((n + absm) * (n + absm - 1) * (n + absm - 2) * Laguerre(n, absm - 3, kSinalSq) + 3 * (n + absm) * (n + absm - 1) * (n + 1) * 
		 Laguerre(n + 1, absm - 3, kSinalSq) + 3 * (n + absm) * (n + 1) * (n + 2) * Laguerre(n + 2, absm - 3, kSinalSq) + 
		 (n + 1) * (n + 2) * (n + 3) * Laguerre(n + 3, absm - 3, kSinalSq) - 2 * (absm - 2) * ((n + absm) * (n + absm - 1) * 
		Laguerre(n, absm - 2, kSinalSq) + 2 * (n + absm) *(n + 1) * Laguerre(n + 1, absm - 2, kSinalSq) + (n + 1) * (n + 2) *
		Laguerre(n + 2, absm - 2, kSinalSq))) / (kSinalSq * kSinal));
}

Complex I2m1(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (std::abs(kSinal) < thresh) {
		if (n % 2 == 0 && (m == 1 || m == -1)) {
			return -I;
		}
		return 0.0;
	}
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;

	return I * (std::exp(- I * absm1 * beta) * I30[n] - std::exp(- I * absm2 * beta) * I31[n]);
}

Complex I31(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			return 2.0 * (n + 1) * m * powm1(n + 1);
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	int sgnm = m >= 0 ? 0 : 1;
	double kSinalSq = kSinal * kSinal;
	return powm1(n + 1 + sgnm) * (kSinal * std::exp(- I * absm1 * beta) * (Laguerre(n, absm + 1, kSinalSq) +
				Laguerre(n - 1, absm + 1, kSinalSq)) + std::exp(- I * absm2 * beta) * 
			((n + absm) * Laguerre(n, absm - 1, kSinalSq) + (n + 1) * Laguerre(n + 1, absm - 1, kSinalSq))/kSinal);
}

Complex I33(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			return 8.0 * (n + 1) * (n + 1) * m * powm1(n + 1);
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	int sgnm = m >= 0 ? 0 : 1;
	double kSinalSq = kSinal * kSinal;
	return powm1(n + sgnm) * (std::exp(- I * absm1 * beta) * (kSinalSq * kSinal * (Laguerre(n, absm + 3, kSinalSq) +
			3.0 * Laguerre(n - 1, absm + 3, kSinalSq) + 3.0 * Laguerre(n - 2, absm + 3, kSinalSq) + Laguerre(n - 3, absm + 3, kSinalSq)) 
			- 2.0 * kSinal * (absm + 2) * (Laguerre(n, absm + 2, kSinalSq) + 2.0 * Laguerre(n - 1, absm + 2, kSinalSq) + 
			Laguerre(n - 2, absm + 2, kSinalSq))) + std::exp(- I * absm2 * beta) * 
		((n + absm) * (n + absm - 1) * (n + absm - 2) * Laguerre(n, absm - 3, kSinalSq) + 3 * (n + absm) * (n + absm - 1) * (n + 1) * 
		 Laguerre(n + 1, absm - 3, kSinalSq) + 3 * (n + absm) * (n + 1) * (n + 2) * Laguerre(n + 2, absm - 3, kSinalSq) + 
		 (n + 1) * (n + 2) * (n + 3) * Laguerre(n + 3, absm - 3, kSinalSq) - 2 * (absm - 2) * ((n + absm) * (n + absm - 1) * 
		Laguerre(n, absm - 2, kSinalSq) + 2 * (n + absm) *(n + 1) * Laguerre(n + 1, absm - 2, kSinalSq) + (n + 1) * (n + 2) *
		Laguerre(n + 2, absm - 2, kSinalSq))) / (kSinalSq * kSinal));
}

Complex I3m1(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (std::abs(kSinal) < thresh) {
		if (n % 2 == 0 && (m == 1 || m == -1)) {
			return -m * 1.0;
		}
		return 0.0;
	}
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	int sgnm = m >= 0 ? 0 : 1;

	return powm1(1 + sgnm) * (std::exp(- I * absm1 * beta) * I30[n] + std::exp(- I * absm2 * beta) * I31[n]);
}

Complex I40(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (n == 0) {
		return 0.0;
	}
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			if (n % 2) {
				return 2.0 * m * (n + 1);
			}
			else {
				return -2.0 * m * n;
			}
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	double sgnm = m >= 0 ? 1 : -1;
	double kSinalSq = kSinal * kSinal;

	return 2.0 * sgnm * (std::exp(- I * absm1 * beta) * powm1(n + 1) * kSinal * Laguerre(n - 1, absm + 1, kSinalSq) +
			std::exp(- I * absm2 * beta) * (1.0 * (n + absm) * I31[n - 1] - 1.0 * n * I31[n]));
}

Complex I42(int n, int m, double kSinal, double beta)
{
	if (n > 0) {
		if (std::abs(kSinal) < thresh) {
			if (m == 1 || m == -1) {
				return 8.0 * n * (n + 1) * m * powm1(n + 1);
			}
			return 0.0;
		}
		int absm = m >= 0 ? m : -m;
		double absm1 = m >= 0 ? m + 1 : m - 1;
		double absm2 =  m >= 0 ? m - 1 : m + 1;
		int sgnm = m >= 0 ? 0 : 1;
		double kSinalSq = kSinal * kSinal;
		return 2.0 * powm1(n + sgnm)* (std::exp(- I * absm1 * beta) * (kSinalSq * kSinal * (Laguerre(n - 1, absm + 3, kSinalSq) +
				2.0 * Laguerre(n - 2, absm + 3, kSinalSq) + Laguerre(n - 3, absm + 3, kSinalSq)) - 2.0 * kSinal * (absm + 2) * 
				   (Laguerre(n - 1, absm + 2, kSinalSq) + Laguerre(n - 2, absm + 2, kSinalSq))) - std::exp(- I * absm2 * beta) * 
			((n + absm - 1) * (n + absm) * Laguerre(n - 1, absm - 1, kSinalSq) + 2 * n * (n + absm) * Laguerre(n, absm - 1, kSinalSq) +
			 n * (n + 1) * Laguerre(n + 1, absm - 1, kSinalSq)) / kSinal);
	}

	return 0.0;
}

Complex I50(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (n == 0) {
		return 0.0;
	}
	if (std::abs(kSinal) < thresh) {
		if (m == 1 || m == -1) {
			if (n % 2) {
				return 2.0 * (n + 1) * I;
			}
			else {
				return -2.0 * n * I;
			}
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 1 : m - 1;
	double absm2 =  m >= 0 ? m - 1 : m + 1;
	double kSinalSq = kSinal * kSinal;

	return 2.0 * I * (std::exp(- I * absm1 * beta) * powm1(n) * kSinal * Laguerre(n - 1, absm + 1, kSinalSq) +
			std::exp(- I * absm2 * beta) * (1.0 * (n + absm) * I31[n - 1] - 1.0 * n * I31[n]));
}

Complex I52(int n, int m, double kSinal, double beta)
{
	if (n > 0) {
		if (std::abs(kSinal) < thresh) {
			if (m == 1 || m == -1) {
				return 8.0 * n * (n + 1) * powm1(n + 1) * I;
			}
			return 0.0;
		}
		int absm = m >= 0 ? m : -m;
		double absm1 = m >= 0 ? m + 1 : m - 1;
		double absm2 =  m >= 0 ? m - 1 : m + 1;
		double kSinalSq = kSinal * kSinal;
		return 2.0 * powm1(n + 1) * I * (std::exp(- I * absm1 * beta) * (kSinalSq * kSinal * (Laguerre(n - 1, absm + 3, kSinalSq) +
				2.0 * Laguerre(n - 2, absm + 3, kSinalSq) + Laguerre(n - 3, absm + 3, kSinalSq)) - 2.0 * kSinal * (absm + 2) * 
				   (Laguerre(n - 1, absm + 2, kSinalSq) + Laguerre(n - 2, absm + 2, kSinalSq))) + std::exp(- I * absm2 * beta) * 
			((n + absm - 1) * (n + absm) * Laguerre(n - 1, absm - 1, kSinalSq) + 2 * n * (n + absm) * Laguerre(n, absm - 1, kSinalSq) +
			 n * (n + 1) * Laguerre(n + 1, absm - 1, kSinalSq)) / kSinal);
	}

	return 0.0;
}

Complex I60(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return -0.25 * m * (1 + (3 + 2 * n) * powm1(n)) * I;
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double sgnm = m >= 0 ? 1.0 : -1.0;
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;

	return I * sgnm * (std::exp(- I * absm1 * beta) * (powm1(n) * Laguerre(n, absm, kSinalSq) - 2.0 * (absm + 1) * I30[n] / kSinal) +
			std::exp(- I * absm2 * beta) * (powm1(n + 1) * Laguerre(n, absm, kSinalSq) + 2.0 * (absm - 1) * I31[n] / kSinal));
}

Complex I62(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return 2.0 * m * (n + 1) * (n + 2) * powm1(n) * I;
		}
		return 0.0;
	}
	
	int absm = m >= 0 ? m : -m;
	double sgnm = m >= 0 ? 1.0 : -1.0;	
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;
	
	return powm1(n) * sgnm * I * (std::exp(- I * absm2 * beta) * ((n + absm) * (n + absm - 1) * Laguerre(n, absm - 2, kSinalSq) +
				2.0 * (n + absm) * (n + 1) * Laguerre(n + 1, absm - 2, kSinalSq) + (n + 1) * (n + 2) * Laguerre(n + 2, absm - 2, kSinalSq)) / kSinalSq  -
				 kSinalSq * std::exp(- I * absm1 * beta) * (Laguerre(n, absm + 2, kSinalSq) + 2.0 * Laguerre(n - 1, absm + 2, kSinalSq) + Laguerre(n - 2, absm + 2, kSinalSq)));
}

Complex I70(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return -0.5 * (1 + (3 + 2 * n) * powm1(n));
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;

	return std::exp(- I * absm1 * beta) * (powm1(n) * Laguerre(n, absm, kSinalSq) - 2.0 * (absm + 1) * I30[n] / kSinal) +
			std::exp(- I * absm2 * beta) * (powm1(n) * Laguerre(n, absm, kSinalSq) - 2.0 * (absm - 1) * I31[n] / kSinal);
}

Complex I72(int n, int m, double kSinal, double beta)
{
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return 4.0 * (n + 1) * (n + 2) * powm1(n + 1);
		}
		return 0.0;
	}
	
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;
	
	return powm1(n + 1) * (std::exp(- I * absm2 * beta) * ((n + absm) * (n + absm - 1) * Laguerre(n, absm - 2, kSinalSq) +
				2.0 * (n + absm) * (n + 1) * Laguerre(n + 1, absm - 2, kSinalSq) + (n + 1) * (n + 2) * Laguerre(n + 2, absm - 2, kSinalSq)) / kSinalSq  +
				 kSinalSq * std::exp(- I * absm1 * beta) * (Laguerre(n, absm + 2, kSinalSq) + 2.0 * Laguerre(n - 1, absm + 2, kSinalSq) + Laguerre(n - 2, absm + 2, kSinalSq)));
}

Complex I451(int n, int m, double kSinal, double beta)
{
	if (n > 0) {
		if (std::abs(kSinal) < thresh) {
			if (m == 0) {
				return 8.0 * n * powm1(n);
			}
			return 0.0;
		}
		int absm = m >= 0 ? m : -m;
		double kSinalSq = kSinal * kSinal;
		return 4.0 * powm1(n) * std::exp(- I * (1.0 * m) * beta) * ((n + absm) * Laguerre(n - 1, absm, kSinalSq) + n * Laguerre(n, absm, kSinalSq));
	}

	return 0.0;
}

Complex I91(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (n == 0) {
		return 0.0;
	}
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return -m * (-1 + (2 * n * n + 4 * n + 1) * powm1(n)) * I;
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double sgnm = m >= 0 ? 1.0 : -1.0;		
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;

	return -2.0 * sgnm * I * (std::exp(- I * absm1 * beta) * powm1(n) * kSinalSq * (Laguerre(n - 1, absm + 2, kSinalSq) +
				Laguerre(n - 2, absm + 2, kSinalSq)) + std::exp(- I * absm2 * beta) * (powm1(n) * ((n + absm) * Laguerre(n - 1, absm, kSinalSq) +
				n * Laguerre(n, absm, kSinalSq)) + 2.0 * (absm - 1) * (1.0 * (n + absm) * I31[n -1] - 1.0 * n * I31[n]) / kSinal));
}

Complex I101(int n, int m, double kSinal, double beta, Complex* I30, Complex* I31)
{
	if (n == 0) {
		return 0.0;
	}
	if (std::abs(kSinal) < thresh) {
		if (m == 2 || m == -2) {
			return -2.0 * (-1 + (2 * n * n + 4 * n + 1) * powm1(n));
		}
		return 0.0;
	}
	int absm = m >= 0 ? m : -m;
	double absm1 = m >= 0 ? m + 2 : m - 2;
	double absm2 =  m >= 0 ? m - 2 : m + 2;
	double kSinalSq = kSinal * kSinal;

	return 2.0 * (std::exp(- I * absm1 * beta) * powm1(n + 1) * kSinalSq * (Laguerre(n - 1, absm + 2, kSinalSq) +
				Laguerre(n - 2, absm + 2, kSinalSq)) + std::exp(- I * absm2 * beta) * (powm1(n) * ((n + absm) * Laguerre(n - 1, absm, kSinalSq) +
				n * Laguerre(n, absm, kSinalSq)) + 2.0 * (absm - 1) * (1.0 * (n + absm) * I31[n -1] - 1.0 * n * I31[n]) / kSinal));
}
