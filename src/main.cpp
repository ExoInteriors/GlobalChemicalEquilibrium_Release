#include "solver.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>

// Function to optimise walkers and return the best parameters and CTOT
int solver::optimize_walkers(int l) {

	for(int i = 0; i < nBestSolutions; ++i){
		bestCTOT[i] = 1.0e10;
	}

	//initialize adam parameters
	eta = initial_eta;
	beta1_t = beta1;
	beta2_t = beta2;

	// Set initial values
	for(int i = 0; i < N * Ndim; ++i) {
		values[i] = initial_values[i];
	}

	setBoundaries(values); 

	// Calculate log of variables
	for(int k = 0; k < N * Ndim; ++k) {
		logvalues[k] = log(values[k]);
	}

	// Calculate the initial pressures
	for(int k = 0; k < N; ++k) {
		P[k] = calculate_P(values + k * Ndim);
	}

	// Calculate initial cost function
	for(int k = 0; k < N; ++k) {
		CTOT[k] = ctot(values + k * Ndim, logvalues + k * Ndim, P[k]);
	}

	//Write initial values
	int er = write_output_to_file(0);
	if(er <= 0) {
		return 0;
	}

	double diffCTOT[N];
	double oldCTOT[N];


	int check_at_iteration = 2000;

	for(int i = 0; i < N; ++i){
		oldCTOT[i] = 1.0e10;
	}

	// Optimization loop
	for(int t = 1; t < nSteps; ++t){
		// Reset SVGD gradients
		for(int k = 0; k < N * Ndim; ++k){
			svgd_grad[k] = 0.0;
		}

		// Compute gradient of the cost function
		for(int k = 0; k < N; ++k){
			grad_ctot(values + k * Ndim, logvalues + k * Ndim, P[k], ctot_grad + k * Ndim);
		}

		// Update gradients
		if(svgd == 1){
			if(useLog == 0){
				SVGD(values);
			}
			else{
			SVGD(logvalues);
			}
		}
		else if(svgd == 2){
			if(useLog == 0){
				SVGD_1(values, N - 1);
			}
			else{
				SVGD_1(logvalues, N - 1);
			}
		}
		else{
			for(int k = 0; k < N * Ndim; ++k){
				svgd_grad[k] = ctot_grad[k];
			}
		}

		// Perform optimization using adamW or adamaxW
		if(method == 1){
			if(useLog == 0){
				adamW_update(values);
			}
			else{
				adamW_update(logvalues);
			}
		} else if(method == 2){
			if(useLog == 0){
				adamaxW_update(values);
			}
			else{
				adamaxW_update(logvalues);
			}
		}


		if(useLog == 0){
			setBoundaries(values);
		}
		else{
			setBoundariesLog(logvalues);
		}


		// Update log and regular variables
		for(int k = 0; k < N * Ndim; ++k){
			if(useLog == 0){
				logvalues[k] = log(values[k]);
			}
			else{
				values[k] = exp(logvalues[k]);
			}
		}

		// Recalculate pressures and cost function
		for(int k = 0; k < N; ++k){
			P[k] = calculate_P(values + k * Ndim);
			CTOT[k] = ctot(values + k * Ndim, logvalues + k * Ndim, P[k]);


			// Update best values
			if(CTOT[k] < bestCTOT[l]){
				bestCTOT[l] = CTOT[k];
				for(int i = 0; i < Ndim; ++i){
					bestValues[l * Ndim + i] = values[k * Ndim + i];
				}
			}


			if(t % check_at_iteration == 0){
				// Take absolute diff between CTOT and oldCTOT
				diffCTOT[k] = abs(oldCTOT[k] - CTOT[k]);
				if(diffCTOT[k] < stopcritera){
					break;
				}
				oldCTOT[k] = CTOT[k];
			}
		}

		// Write to file every outputInterval
		er = write_output_to_file(t);
		if(er <= 0){
			return 0;
		}

	}//end of optimization loop

	// Store the best solutions from the last iteration
	std::vector<std::pair<double, std::vector<double>>> bestSolutions; // Pair of CTOT and corresponding values

	for(int k = 0; k < N; ++k){
		// Store CTOT and corresponding values as a pair
		std::vector<double> solution(values + k * Ndim, values + (k + 1) * Ndim); // Extract values[k]
		bestSolutions.emplace_back(CTOT[k], solution);
	}

	// Sort by CTOT (ascending order) to get the lowest values
	std::sort(bestSolutions.begin(), bestSolutions.end(),
	[](const std::pair<double, std::vector<double>> &a, const std::pair<double, std::vector<double>> &b) {
		return a.first < b.first;
	});

	// Store the top nBestSolutions in bestCTOT and bestValues
	for(int i = 0; i < nBestSolutions; ++i){
		bestCTOT[i] = bestSolutions[i].first;
		for(int d = 0; d < Ndim; ++d){
			bestValues[i * Ndim + d] = bestSolutions[i].second[d];
		}
	}

	return 1;
}

//int main(int argc, char *argv[]) {
int main() {

	// Initilize solver
	solver S;
	int er = S.readParam();

	if(er <= 0) {
		return 0;
	}

	//Allocate memory
	S.allocate();


	// Read chemical input
	er = S.readcheminput();
	if(er <= 0) {
		return 0;
	}

	// Read Gibbs Energies
	er = S.readGibbs();
	if(er <= 0) {
		return 0;
	}

	//Set SVGD bandwidth
	for(int i = 0; i < 26; ++i){
		S.h[i] = 0.5;
	}
	for(int i = 26; i < S.Ndim; ++i){
		S.h[i] = 0.5 * 10000.0;
	}


	// Generate random initial values
	std::default_random_engine generator;
	generator.seed(std::time(0));
	//generator.seed(0);

	if(S.useICFile == 0){

		for(int k = 0; k < S.N; ++k){
			for(int i = 0; i < S.Ndim; ++i){
				std::uniform_real_distribution<double> distribution(S.prior_min[i], S.prior_max[i]);
				S.initial_values[k * S.Ndim + i] = distribution(generator);
			}
		}
	}
	else{
		//Read initial values for a single walker from file
		printf("Read initial conditions from file %s\n", S.inputFilename);
		er = S.readIC();
		if(er <= 0){
			return 0;
		}
		//Add perturbation to other walkers
		for(int k = 1; k < S.N; ++k){
			for(int i = 0; i < S.Ndim; ++i){
				std::normal_distribution<double> distribution(0.0, 0.001 * (S.prior_max[i] - S.prior_min[i]));

				S.initial_values[k * S.Ndim + i] = S.initial_values[i] + distribution(generator);
			}
		}
	}

	S.setBoundaries(S.initial_values);

	if(S.svgd == 2) {
		S.N = 1;
	}

	// Vector to store the top solutions
	std::vector<std::pair<double, std::vector<double>>> bestSolutions;


	//Start iteration loop
	for(int iter = 0; iter < S.nIterations; ++iter){
		printf("%d ", iter+1);

		// Container for all solutions from the current iteration
		std::vector<std::pair<double, std::vector<double>>> solutionPool;

		// Run optimization for each of the nBestSolutions
		for(int l = 0; l < (iter == 0 ? 1 : S.nBestSolutions); ++l){
			if(iter > 0){
				// Subsequent iterations: perturb the l-th best solution
				for(int k = 0; k < S.N; ++k){
					for (int i = 0; i < S.Ndim; ++i) {
						std::normal_distribution<double> distribution(0.0, 1.0e-1 * bestSolutions[l].second[i]);
						double random = distribution(generator);
						S.initial_values[k * S.Ndim + i] = bestSolutions[l].second[i] + random;
					}
				}
			}

			// Optimize the walkers
			int er = S.optimize_walkers(l);
			if (er <= 0) return 0;

			// Collect solutions
			for(int i = 0; i < S.nBestSolutions; ++i){
				std::vector<double> solution(S.bestValues + i * S.Ndim, S.bestValues + (i + 1) * S.Ndim);
				solutionPool.emplace_back(S.bestCTOT[i], solution);
			}
		}

		// Sort all collected solutions by CTOT
		std::sort(solutionPool.begin(), solutionPool.end(),
		[](const std::pair<double, std::vector<double>> &a, const std::pair<double, std::vector<double>> &b) {
		return a.first < b.first;
		});

		// Keep only the top nBestSolutions
		bestSolutions.clear();
		for(int i = 0; i < S.nBestSolutions; ++i){
			bestSolutions.push_back(solutionPool[i]);
		}

		// Print the best solutions for this iteration
		for(int i = 0; i < S.nBestSolutions; ++i){
			printf("%g ", bestSolutions[i].first);
			// Print parameters
			for(int j = 0; j < S.Ndim; ++j){
				printf("%g ", bestSolutions[i].second[j]);
			}
			printf("\n");
		}

		//needs to be updated
		if(S.svgd == 2) {
			for(int i = 0; i <= S.Ndim-3; ++i){
				std::normal_distribution<double> distribution(0.0, 1.0e-1 * S.values[(S.N - 1) * S.Ndim + i]);
				double random = distribution(generator);
				S.initial_values[S.N * S.Ndim + i] = S.values[(S.N - 1) * S.Ndim + i] + random;
			}
			for(int i = S.Ndim-2; i <= S.Ndim; ++i){
				std::normal_distribution<double> distribution(0.0, 1.0e-1 * S.values[(S.N - 1) * S.Ndim + i]);
				double random = distribution(generator);
				S.initial_values[S.N * S.Ndim + i] = S.values[(S.N - 1) * S.Ndim + i] + random;
			}
			++S.N;
		}
	}//end iteration loop

	return 0;
}
