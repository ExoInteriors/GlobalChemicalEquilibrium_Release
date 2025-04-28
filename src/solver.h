#ifndef SOLVER_H
#define SOLVER_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <string>


#include <ctime>	// Required for seeding the random number generator

// *****************************************************************************************
// Author: Simon Grimm
// Date: February 2025
// *****************************************************************************************


class solver{


public:


	// *****************************************************************************
	//Chemical parameters
	#include "../Standard_Version/parameters.h"
	//The definition of the Parameters array, including "nSi, nMg, nO, ...."
	//is moved to the automatically generated file "parameters,h" and included here
	// *****************************************************************************



	double *P;
	double *GRT_T;
	double *prior_min;		//minimal value for initial sampling
	double *prior_max;		//maximal value for initial sampling
	double *bound_min;		//minimal value
	double *bound_max;		//maximal value
	double *logbound_min;		//minimal value
	double *logbound_max;		//maximal value

	double *CTOT;			// cost function
	double *ctot_grad;		// gradients of the cost function
	// ****************************************************************************

	int N = 1;			// Number of walkers
	int outputInterval = 1;
	int nSteps = 1;

	int nIterations = 1;
	int nBestSolutions = 1;


	// SVGD parameters
	double *h;			// bandwidth
	double *svgd_grad;		// SVGD gradient for each parameter
	double f1 = 1.0;		// Attractive force
	double f2 = 1.0;		// Repulsive force
	int svgd = 0;			// Use SVGD or not

	double *bestCTOT;


	// AdamW parameters
	double eta = 0.4;		// Learning rate
	double initial_eta = 0.4;	// initial Learning rate
	double eps = 1.0e-6;		// Small constant for numerical stability
	double beta1 = 0.9;
	double beta2 = 0.999;
	double beta1_t = beta1;
	double beta2_t = beta2;
	double lambda_reg = 0.0001;	// Regularization parameter for AdamW

	double stopcritera = 1.0e-15;

	// Optimizer
	int method = 2;			// 1: AdamW, 2: AdamaxW

	bool reset_best_values;

	int useICFile = 0;
	char inputFilename[160];
	char outputFilename[160];



	double *m;
	double *v;
	double *d;

	double *initial_values;
	double *bestValues;
	double *values;
	double *logvalues;



	void allocate();
	int readParam();
	int readcheminput();
	int readIC();
	void setBoundaries(double *);
	void setBoundariesLog(double *);
	void adamW_update(double *);
	void adamaxW_update(double *);
	void SVGD(double *);
	void SVGD_1(double *, int);
	double calculate_P(double *);
	double ctot(double *, double *, double);
	void grad_ctot(double *, double *, double, double *);
	int optimize_walkers(int);

	int write_output_to_file(int);


};


#endif //SOLVER_H

