#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <format>

#if __cplusplus < 202002L
#error "This file requires C++20 or later."
#endif


int main(){
	std::string line;
	std::ifstream inFile("data_unix.txt"); // Input file
	std::ofstream outFile("output.txt"); // Output file

	// Check if input file exists
	if (inFile.is_open()){
		int exponent = 0;
		double num = 0.;
		// Iterate through all lines of input file
		while(getline(inFile,line)){
			// Simple method using built-in std::scientific
//			std::cout << std::scientific << std::stod(line) << std::endl;
//			outFile << std::scientific << std::stod(line) << ", ";

			// Manual method
			// Calculate exponent
			exponent = line.substr(0, line.find(".")).size() - 1;
			// Convert str to double and divide by 1e exponent
			num = std::stod(line);
			num /= std::pow(10, exponent);
			// Format number to scientific notation and print to cli as well as to outFile
//			std::cout << std::format("{:.6f}", num) << "e" << std::format("{:02}", exponent) << std::endl;
			outFile << std::format("{:.6f}", num) << "e" << std::format("{:02}", exponent) << ",";
		}
		// Close in- and output files
		std::cout << "Output written to output.txt." << std::endl;
		inFile.close();
		outFile.close();
	}

	else std::cout << "Unable to open file";

	return 0;
}



