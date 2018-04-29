
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

double	h	= 1;
double	m 	= 1;
double	w	= 1;
int		D 	= 8;
int 	N 	= 10;
double	XS	= 0.75;	
double	X0	= -4;
double	XD	= 4;
double 	PI  = 3.14159265358;	
double	T0 	= 2 * PI;
double	T 	= T0/16;
double	DEL_T 	= T0 / 128;
double	DEL_X	= (XD - X0)/D;
double	ALPHA	= 2;

double sin(double x);
double cos(double x);
double exp(double x);


// ===================================
// DECLARE FUNCTIONS
// ===================================

matrix returnPropagator(int n);
vector<cdouble> returnWaveFunction(int n);
double avgPosition(vector<cdouble> wf);

cdouble printMatrix(matrix k);
vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);


// ===================================
// RUN QUESTIONS 1-5
// ===================================


int main() 
{
	// Question 1
	matrix K = returnPropagator(N);
	// printMatrix(K);

	// Question 2
	vector<cdouble> wf = returnWaveFunction(N);
	// printMatrix(wf);

	// Question 3


	// Question 4
	double avg_pos = avgPosition(wf);
	cout<<avg_pos<<endl;

	// vector<cdouble> wf_iter;
	// double avg_pos_iter;
	// for (int i=0; i<N; i++) {
	// 	wf_iter = returnWaveFunction(i);
	// 	avg_pos_iter = avgPosition(wf_iter);
	// 	cout<<avg_pos_iter<<endl;
	// }

	return 0;
}


// ===================================
// QUESTION 1 FUNCTIONS
// ===================================

double x(int d)
{
	return X0 + d*DEL_X;
}


vector<cdouble> phiInit()
{
	vector<cdouble> phi0(D,0);

	for (int i=0; i<D; i++) {
		phi0[i] = pow(ALPHA/PI, .25) * exp(-ALPHA * pow((x(i) - XS), 2) /2);
	}

	return phi0;
}


matrix returnPropagator(int n) 
{
	vector<cdouble> Kep(D*D,0);

	double A = 1; 
	double V = .5*m*pow(w,2);
	cdouble I(0,1);
	cdouble efactor;

	// Generate elementary Kep matrix
	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {

			efactor = I*DEL_T/h * (.5*m*pow(x(j)-x(i),2) / pow(DEL_T,2) - V*(x(j)-x(i))/2);
			// cout<<exp(efactor)<<endl;

			Kep[i*D + j] = exp(efactor);
		}
	}

	vector<cdouble> K(D*D,0);
	K = Kep;

	// Multiply to the nth power to propagate by n time steps
	for (int p=0; p<n; p++) {
		K = matrixMatrixMultiply(K, Kep);
	}

	// Multiply by step size factors
	double factor = pow(DEL_X, n-1) * DEL_T;

	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {
			K[i*D + j] *= factor;
		}
	}

	return K;
}


// ===================================
// QUESTION 2 FUNCTIONS
// ===================================

vector<cdouble> normalizeWF(vector<cdouble> wf)
{

	vector<cdouble> norm_wf(D,0);
	cdouble sum;

	// Square wave function: multiply by complex conjugate
	for (int i; i<D; i++) {
		sum += wf[i] * conj(wf[i]);
	}

	// Divide by normalization factor such that sum(phi*phi) = 1
	for (int i=0; i<D; i++) {
		norm_wf[i] = wf[i] / sum;
	}

	return norm_wf;
}


vector<cdouble> returnWaveFunction(int n)
{
	// Get discretized initial wave function array
	vector<cdouble> phi0 = phiInit();

	// Get propagator matrix after n timesteps
	vector<cdouble> Kn = returnPropagator(n);

	// Multiply propagator matrix and initial wave function
	vector<cdouble> wf = matrixVectorMultiply(Kn, phi0);

	// Normalize wave function
	vector<cdouble> norm_wf = normalizeWF(wf);

	return norm_wf;
}


// ===================================
// QUESTION 3 FUNCTIONS
// ===================================

cdouble avgEnergy()
{
	return 0;
}


cdouble avgKinetic()
{
	return 0;
}


cdouble avgPotential()
{
	return 0;
}


// ===================================
// QUESTION 4 FUNCTIONS
// ===================================

double avgPosition(vector<cdouble> wf)
{
	int n = wf.size();
	cdouble sum;

	for (int i=0; i<n; i++) {
		sum += wf[i];
	}

	double avg = sum.real() / n;

	return avg;
}


// ===================================
// MATRIX/VECTOR OPERATION FUNCTIONS
// ===================================

cdouble printMatrix(matrix k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}

vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    //Create a length D^2 array of doubles, filled with the value 0.
    vector<cdouble> matrixnew(D*D,0);

    //cblas zgemm takes in three matrices: A,B,C. It stores in the value C
    //the matrix alpha*AB+beta*C. In this case beta=0 and alpha=1, so we just
    //do matrix multiplication AB.

    cdouble mult(1, 0.0);
    cdouble czero(0.0, 0.0);

    //use complex matrix matrix multiply.
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,D,D,D, &mult, &m1[0], D, &m2[0], D, &czero, &matrixnew[0], D);
    
    //return the vector object.
    return matrixnew;
}


vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

	vector<cdouble> prod(D,0);

    for (int i=0; i<D; i++) {
    	for (int j=0; j<D; j++) {
    		prod[i] += m1[i*D + j] + m2[j];
    	}
    }

    return prod;
}


vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

	vector<cdouble> prod(D,0);

    for (int i=0; i<D; i++) {
    	prod[i] += m1[i] + m2[i];
    }

    return prod;
}

