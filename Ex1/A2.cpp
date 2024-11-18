#include <iostream>
#include <cmath>

#if __cplusplus < 202002L
#error "This file requires C++20 or later."
#endif


/**
 * Calculates the length of the given vector
 *
 * @param v Three dimensional vector to calculate length of
 * @return Length of vector
 */
double vectorLength(const double* v) {
	return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


/**
 * Normalizes given vector (in-place)
 *
 * @param v Three dimensional vector to normalize
 */
void normalizeVector(double* v) {
	double length = vectorLength(v);
	if (length > 0) {
		v[0] /= length;
		v[1] /= length;
		v[2] /= length;
	}
}

/**
 * Checks if the given vectors are linear dependent in R3
 *
 * @param a 3-dim vector as double array
 * @param b 3-dim vector as double array
 * @param c 3-dim vector as double array
 * @return Boolean if vectors are linear dependent
 */
bool LinearDependent(const double* a, const double* b, const double* c) {
	// Vectors are linear independent if space spanned by those vectors has a orthogonal base
	// This means that Matrix (v1, v2, v3) needs to be full rank -> regular
	// -> det(v1, v2, v3) = 0 and in R3 det(v1, v2, v3) = v1 * (v2 x v3) = 0

	// Working with numerical values -> normed vectors can achieve better accuracy without doing a SVD

	// Create copies of the vectors to normalize
	double a_norm[3] = {a[0], a[1], a[2]};
	double b_norm[3] = {b[0], b[1], b[2]};
	double c_norm[3] = {c[0], c[1], c[2]};

	// Normalize the vectors
	normalizeVector(a_norm);
	normalizeVector(b_norm);
	normalizeVector(c_norm);

	// Calculate the determinant of the matrix of the normalized vects
	double det = a_norm[0] * (b_norm[1] * c_norm[2] - b_norm[2] * c_norm[1]) -
	         a_norm[1] * (b_norm[0] * c_norm[2] - b_norm[2] * c_norm[0]) +
	         a_norm[2] * (b_norm[0] * c_norm[1] - b_norm[1] * c_norm[0]);

	// Check toleracne criteria
	return std::abs(det) < 1e-9;
}

int main(){
	// Test 1
	const double a1[3] = {1, 2, 3};
	const double b1[3] = {-.6, 3., 5};
	const double c1[3] = {1, 2.2, -2.4};

	// Test 2
	const double a2[3] = {1./3, 2./3, 1};
	const double b2[3] = {-1e8, 5e8, 0};
	const double c2[3] = {0, 1, 3./7};

	//Test 3
	const double a3[3] = {1./3, 2./3, 1};
	const double b3[3] = {-1e-8, 7e-8, 0};
	const double c3[3] = {0, 1, 3./7};

	// Test cases and print to cli
	std::cout << "Test 1 (independent): " << (LinearDependent(a1, b1, c1) ? "dependent" : "independent") << std::endl;
	std::cout << "Test 2 (dependent):   " << (LinearDependent(a2, b2, c2) ? "dependent" : "independent") << std::endl;
	std::cout << "Test 3 (independent): " << (LinearDependent(a3, b3, c3) ? "dependent" : "independent") << std::endl;
	return 0;
}
