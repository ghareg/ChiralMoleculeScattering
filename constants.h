#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double Pi = 3.14159;

const int nmax = 2;
const int mmax = 2;

const double Lz = 15;
const double a = 3;

const double mu1 = 1.2;
const double mu2 = -1.9;
const double mu3 = 0.36;
const double DipConst = 4 * 3.33 * 1.6 * 9.1 / (Pi * 8.85 * 1.05 * 1.05);


const double F1 = DipConst * mu1 / (Lz * Lz * Lz);
const double F2 = DipConst * mu2 / (2 * 2 * 2);
const double F3 = DipConst * mu3 / (2 * 2 * 2);

const double EnMult = 2 * 9.1 * 1.6 * a * a * 1E-5 / (1.05 * 1.05);

const double alpha0 = 1.05 * 1.05 * 0.01 / (4 * 9.1 * 9.1 * 9 * a * a); 
const double alphaSOC = alpha0 * 1E6; //alpha / (a * a); alpha in A^2 

const int ND = 100;
const double EMax = 1000;

const int nCore = 4;
const double thresh = 1E-6;


#endif
