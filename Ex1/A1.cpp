#include <iostream>
#include <cmath>

#include <format>
#include <iomanip>	// for io manipulation i.e. setw

#if __cplusplus < 202002L
#error "This file requires C++20 or later."
#endif

/**
 * Computes the taylor series of log(x) with x0 = 1 and N steps
 *
 * @param x Value for which log is approximated for
 * @param N Final value to calculate sum to
 * @return Approximated log(x) using taylor series
 */
double taylor_log(double x, int N){
	double sum = 0;
	for(int n = 1; n <= N; n++){
		sum += pow(-1, n+1) * pow(x-1, n) / n;
	}
	return sum;
}


/**
 * Computes the series1 of log(x) with N steps
 *
 * @param x Value for which log is approximated for
 * @param N Final value to calculate sum to
 * @return Approximated log(x) using series one
 */
double series1_log(double x, int N){
	double sum = 0;
	for(int n = 1; n <= N; n++){
		sum += pow((x-1) / x, n) / n;
	}
	return sum;
}

/**
 * Computes the series1 of log(x) with N steps
 *
 * @param x Value for which log is approximated for
 * @param N Final value to calculate sum to
 * @return Approximated log(x) using series two
 */
double series2_log(double x, int N){
	double sum = 0;
	for(int n = 1; n <= N; n++){
		sum += pow((x-1) / (x+1), 2*n-1) / (2*n - 1);
	}
	return 2*sum;
}


/**
 * Iterates over N until the approximated log is in tolerance
 *
 * @param x Value for which log is approximated for
 * @param log_func Function to use for approximation
 * @param true_log Real log value using log(x)
 * @param tol Toleracnce for approximation criteria
 * @return Number of finale value needed to be in tolerance for approximation
 */
int find_N(double x, double (*log_func)(double, int), double true_log, double tol){
	int N = 1;
	while(fabs(log_func(x, N) - true_log) > tol){
		N++;
	}
	return N;
}

int main() {
	double x_vals[] = {1.1, 0.51, .1, 1.9};
	double tol = 1e-6;

	// Iterate through x values
	for(double x : x_vals){
		double true_log = std::log(x);

		// Format header of table
		std::cout << std::string(34, '=') << std::endl;
		std::cout << std::left << "|" << std::setw(8) << std::format("x = {:.2f}", x) << std::right<< std::setw(24) <<std::format("log(x) = {:.7f}", true_log) << "|" << std::endl;
		std::cout << "|" << std::string(32, '-') << "|" << std::endl;
		std::cout << std::left << "|" << std::setw(10) << "Series" << " | " << std::setw(5) << "N" << " | " << std::setw(10) << "series(x)" << " |" << std::endl;
		std::cout << "|" << std::string(32, '=') << "|" <<std::endl;

		// If x is in convergence interval use series to approximate log(x)
		if (x > 0 && x <= 2){
			int N_taylor = find_N(x, taylor_log, true_log, tol);
			std::cout << std::left << "|" << std::setw(10) << "Taylor" << " | " << std::setw(5) << N_taylor << " | " << std::setw(10) << std::format("{:.7f}", taylor_log(x, N_taylor)) << " |" << std::endl;
		}
		if (x > .5){
			int N_series1 = find_N(x, series1_log, true_log, tol);
			std::cout << std::left << "|" << std::setw(10) << "Series1" << " | " << std::setw(5) << N_series1 << " | " << std::setw(10) << std::format("{:.7f}", series1_log(x, N_series1)) << " |" << std::endl;
		}
		if (x > 0){
			int N_series2 = find_N(x, series2_log, true_log, tol);
			std::cout << std::left << "|" << std::setw(10) << "Series2" << " | " << std::setw(5) << N_series2 << " | " << std::setw(10) << std::format("{:.7f}", series2_log(x, N_series2)) << " |" << std::endl;
		}
		std::cout << std::string(34, '=') << std::endl << std::endl;
	}
	return 0;
}
