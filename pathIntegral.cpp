// ===================================
// JESSICA BIRKY (A13002163)
// ===================================

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

int		D 	= 601; 				// number of grid points
int 	N 	= 8; 				// number of time steps
double	h	= 1;
double	m 	= 1;
double	w	= 1;
double	XS	= 0.75;				// center of initial condition gaussian
double	X0	= -4;				// min x range
double	XD	= 4;				// max s range
double 	PI  = 3.14159265358;	
double	T0 	= 2 * PI;
int 	P 	= 16;
double	T 	= T0/P;				// period
double	DEL_T 	= T0 / 128;		// time step size
double	DEL_X	= (XD - X0)/D;	// spatial step size
double	ALPHA	= 2;

double sin(double x);
double cos(double x);
double exp(double x);


// ===================================
// DECLARE FUNCTIONS
// ===================================

vector<cdouble> phiInit();
matrix 			returnKep();
matrix 			returnPropagator(int n);					// return propagator matrix at time step n
vector<cdouble> returnWaveFunction(int n);					// return wave function at time step n
void			saveWaveFunctions(int n);					// save wave functions, probability amplitudes, and expected values at each time step
vector<double> 	returnProbability(vector<cdouble> wf);		// return wave function squared
vector<cdouble> normalizeWF(vector<cdouble> wf);			// normalize wave function such that sum(wf* wf) = 1

double 	avgPosition(vector<cdouble> wf); 					// return <x> for each time step up to n
double 	avgKinetic(vector<cdouble> wf);	 					// return <K> for each time step up to n
double 	avgPotential(vector<cdouble> wf); 					// return <V> for each time step up to n

cdouble printMatrixC(matrix k); 							// print complex matrix
double  printMatrixR(vector<double> k); 					// print real matrix
void 	saveFileC(vector<cdouble> vec, string save_name);	// save complex matrix
void 	saveFileR(vector<double> vec,  string save_name);	// save real matrix

vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);
vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2);


// ===================================
// RUN QUESTIONS 1-5
// ===================================


int main() 
{
	clock_t begin = clock();

	// ========================================= 

	// K = returnPropagator(N);
	// printMatrixC(K);


	// vector<cdouble> wf = returnWaveFunction(N);
	// printMatrixC(wf);

	saveWaveFunctions(N);

	// ========================================= 

	// matrix K;
	// string mat_name;
	// for (int i=0; i<N; i++) {
	// 	K = returnPropagator(i);
	// 	mat_name = "prop_matrix/Kprop" + to_string(i) + ".dat";
	// 	saveFileC(K, mat_name);
	// }

	// vector<cdouble> phi0 = phiInit();
	// string phi0_name = "wave_func/phi0.dat";
	// saveFileC(phi0, phi0_name);

	// vector<double> phi_sq0 = returnProbability(phi0); 
	// string phi0_sq = "wave_prob/phi_sq0.dat";
	// saveFileR(phi_sq0, phi0_sq);

	// matrix K1 = returnPropagator(1);
	// matrix phi1 = matrixVectorMultiply(K1, phi0);
	// phi1 = normalizeWF(phi1);
	// vector<double> phi_sq1 = returnProbability(phi1); 
	// string phi1_sq = "wave_prob/phi_sq1.dat";
	// saveFileR(phi_sq1, phi1_sq);
	
	// =========================================

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / pow(10,6);
	cout<<"Time elapsed: "<<elapsed_secs<<" sec"<<endl;

	return 0;
}


// ===================================
// CREATE PROPAGATOR MATRIX
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


matrix returnKep()
{
	matrix Kep(D*D,0);

	// cdouble A(sqrt(PI * h * DEL_T), sqrt(PI * h * DEL_T)); 
	double A = 1; // normalize wave function later
	double V = .5*m*pow(w,2);
	cdouble I(0,1);
	cdouble efactor;

	// Generate elementary Kep matrix
	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {

			efactor = I*DEL_T/h * (.5*m*pow(x(j)-x(i),2) / pow(DEL_T,2) - V*(x(j)-x(i))/2);
			// cout<<exp(efactor)<<endl;

			Kep[i*D + j] = exp(efactor) / A;
		}
	}

	return Kep;
}


matrix returnPropagator(int n) 
{
	matrix Kep = returnKep();
	matrix K(D*D,0);
	// matrix K = Kep;

	// Set K equal to identitity matrix for time step 0
	for (int i=0; i<D; i++) {
		for (int j=0; j<D; j++) {
			if (i==j) {
				K[i*D + j] = 1;
			} else {
				K[i*D + j] = 0;
			}
		}
	}

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
// COMPUTE WAVE FUNCTION VECTOR
// ===================================

vector<cdouble> normalizeWF(vector<cdouble> wf)
{

	vector<cdouble> norm_wf(D,0);
	double sum = 0;

	// Square wave function: multiply by complex conjugate
	for (int i=0; i<D; i++) {
		sum += (wf[i] * conj(wf[i])).real() * DEL_X;
	}

	// Divide by normalization factor such that sum(phi*phi) = 1
	for (int i=0; i<D; i++) {
		norm_wf[i] = wf[i] / sqrt(sum);
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


vector<double> returnProbability(vector<cdouble> wf) 
{
	vector<double> wf_sq(D,0);

	for (int i=0; i<D; i++) {
		wf_sq[i] = (conj(wf[i]) * wf[i]).real();
	}

	return wf_sq;
}


void saveWaveFunctions(int n)
{
	vector<cdouble> wf;
	vector<double> wf_sq;
	string sname0;
	string sname1;
	vector<double> avg_pos(n,0);
	vector<double> avg_eng(n,0);
	vector<double> avg_kin(n,0);
	vector<double> avg_pot(n,0);

	for (int i=0; i<n; i++) {

		// Save wave function vectors
		wf = returnWaveFunction(i);
		sname0 = "wave_func/phi" + to_string(i) + ".dat";
		saveFileC(wf, sname0);

		// Save probability amplitudes
		wf_sq = returnProbability(wf);
		sname1 = "wave_prob/phi_sq" + to_string(i) + ".dat";
		saveFileR(wf_sq, sname1);

		// Compute average values
		avg_pos[i] = avgPosition(wf);
		avg_pot[i] = avgPotential(wf);
		avg_kin[i] = avgKinetic(wf);
		avg_eng[i] = avg_pot[i] + avg_kin[i];
	}

	// Save average values
	string sname2 = "expected/avg_pos.dat";
	string sname3 = "expected/avg_pot.dat";
	string sname4 = "expected/avg_kin.dat";
	string sname5 = "expected/avg_eng.dat";
	saveFileR(avg_pos, sname2);
	saveFileR(avg_pot, sname3);
	saveFileR(avg_kin, sname4);
	saveFileR(avg_eng, sname5);
}


// ===================================
// EXPECTED VALUES
// ===================================

double avgPosition(vector<cdouble> wf)
{
	double avg_pos;
	
	for (int i=0; i<D; i++) {
		avg_pos += (conj(wf[i]) * x(i) * wf[i]).real();
	}

	return avg_pos;
}


double avgKinetic(vector<cdouble> wf)
{
	double avg_kin;

	for (int i=0; i<D; i++) {
		if (i == 0) {
			avg_kin += avg_kin += (-h*h / (2*m)) * ((conj(wf[i]) * wf[i+1] - 2.*wf[i])).real();
		} else if (i == D) {
			avg_kin += (-h*h / (2*m)) * ((conj(wf[i]) - 2.*wf[i] + wf[i-1])).real();
		} else {
			avg_kin += (-h*h / (2*m)) * ((conj(wf[i]) * wf[i+1] - 2.*wf[i] + wf[i-1])).real();
		}
	}

	return avg_kin;
}


double avgPotential(vector<cdouble> wf)
{
	double avg_pot;

	for (int i=0; i<D; i++) {
		avg_pot = (conj(wf[i]) * .5*m*w*pow(x(i),2) * wf[i]).real();
	}

	return avg_pot;
}


// ===================================
// MATRIX/VECTOR OPERATION FUNCTIONS
// ===================================

cdouble printMatrixC(matrix k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


double printMatrixR(vector<double> k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


void saveFileC(vector<cdouble> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
}


void saveFileR(vector<double> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
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
