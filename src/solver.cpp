#include "solver.h"



//Allocate all arrays
//N is the number of chains,
//Ndim is the number of equations
void solver::allocate(){

	GRT_T = (double*)malloc(nGRT * sizeof(double));
	prior_min = (double*)malloc(Ndim * sizeof(double));
	prior_max = (double*)malloc(Ndim * sizeof(double));
	bound_min = (double*)malloc(Ndim * sizeof(double));
	bound_max = (double*)malloc(Ndim * sizeof(double));
	logbound_min = (double*)malloc(Ndim * sizeof(double));
	logbound_max = (double*)malloc(Ndim * sizeof(double));

	for(int i = 0; i < Ndim; ++i){
		prior_min[i] = 0.0;
		prior_max[i] = 0.0;
		bound_min[i] = 0.0;
		bound_max[i] = 0.0;
	}
	
	P = (double*)malloc(N * sizeof(double));
	CTOT = (double*)malloc(N * sizeof(double));
	
	h = (double*)malloc(Ndim * sizeof(double));
	initial_values = (double*)malloc(N * Ndim * sizeof(double));
	bestValues = (double*)malloc(nBestSolutions * Ndim * sizeof(double));

	bestCTOT = (double*)malloc(nBestSolutions * sizeof(double));
	svgd_grad = (double*)malloc(Ndim * N * sizeof(double));
	m = (double*)malloc(Ndim * N * sizeof(double));
	v = (double*)malloc(Ndim * N * sizeof(double));
	d = (double*)malloc(Ndim * N * sizeof(double));
	values = (double*)malloc(Ndim * N * sizeof(double));
	logvalues = (double*)malloc(Ndim * N * sizeof(double));
	ctot_grad = (double*)malloc(Ndim * N * sizeof(double));

	for(int k = 0; k < N; ++k){
		for(int i = 0; i < Ndim; ++i){
			m[k * Ndim + i] = 0.0;
			v[k * Ndim + i] = 0.0;
		}
	}


}

// Function to write output to a file
int solver::write_output_to_file(int t) {

	if (t % outputInterval == 0 || t == nSteps) {
		FILE *file = fopen(outputFilename, "a"); // Open file in append mode
		if (file == NULL) {
			perror("Error opening file");
			return 0;
		}

		// Write output if the condition is met
		for (int k = 0; k < N; ++k) {
		//for (int k = N -1; k < N; ++k) {
			fprintf(file, "%d %d %g ", t, k, CTOT[k]);
			for (int i = 0; i < Ndim; ++i) {
				fprintf(file, "%g ", values[k * Ndim + i]);
			}
			fprintf(file, "\n");
		}
		fclose(file); // Close the file
	}
	return 1;
}

int solver::readcheminput(){

	printf("Start reading chemical inputfile\n");

	FILE *paramfile;
	paramfile = fopen("chem_input.dat", "r");

	if(paramfile == NULL){
		printf("Error, chem_input.dat file does not exist\n");
		return 0;
	}

	//------------------------------------------------

	char sp[160];
	char st[160];


	char skip[8];
	int er; 
	int pCount = 0;
	int gCount = 0;

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the param.dat file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(paramfile);
			if(c == EOF){
				break;
			}
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
			if(c == '\n'){
				//blank line
				i = -1;
				continue;
			}
		}
		if(c == EOF) break;
//printf("j = %d %s\n", j, sp);

		//Read GRT array
		int fg = 0;
		for(int p = 0; p < nGRT; ++p){
			sprintf(st, "GRT_%d =", p);

			if(strcmp(sp, st) == 0){
				er = fscanf (paramfile, "%lf", &GRT_T[p]);
				if(er <= 0){
					printf("Error: GRT_%d is not valid!\n", p);
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fg = 1;
				++gCount;
			}
		}
		if(fg == 1) continue;


		//Read Parameters array
		int fp = 0;
		for(int p = 0; p < nParameters; ++p){
			sprintf(st, "%s =", ParameterNames[p].c_str());
			if(strcmp(sp, st) == 0){
				er = fscanf (paramfile, "%lf", &Parameters[p]);
//printf("%s %g\n", st, Parameters[p]);
				if(er <= 0){
					printf("Error: %s is not valid!\n", ParameterNames[p].c_str());
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fp = 1;
				++pCount;
			}
		}
		if(fp == 1) continue;

		//Read Bound array
		int fbound = 0;
		for(int p = 0; p < Ndim; ++p){
			sprintf(st, "bound_%s =", VariableNames[p].c_str());
			if(strcmp(sp, st) == 0){
				er = fscanf (paramfile, "%lf %s %lf", &bound_min[p], skip, &bound_max[p]);
//printf("%s %g %g| %d\n", st, bound_min[p], bound_max[p], p);
				logbound_min[p] = log(bound_min[p]);
				logbound_max[p] = log(bound_max[p]);
				if(er <= 0){
					printf("Error: bound_%d is not valid!\n", p);
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fbound = 1;
			}
		}
		if(fbound == 1) continue;

		//Read Priors array
		int fprior = 0;
		for(int p = 0; p < Ndim; ++p){
			sprintf(st, "prior_%s =", VariableNames[p].c_str());
			if(strcmp(sp, st) == 0){
				er = fscanf (paramfile, "%lf %s %lf", &prior_min[p], skip, &prior_max[p]);
//printf("%s %g %g| %d\n", st, prior_min[p], prior_max[p], p);
				if(er <= 0){
					printf("Error: prior_%d is not valid!\n", p);
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fprior = 1;
			}
		}
		if(fprior == 1) continue;

		
		printf("Error: chem_input.dat file is not valid! %s\n", sp);
		return 0;
	}
	fclose(paramfile);
	if(pCount != nParameters){
		printf("Error, number of parameters in chem_input.dat file does not agree with nParameters\n");	
		return 0;
	}
	if(gCount != nGRT){
		printf("Error, number of GRT terms in chem_input.dat file does not agree with nGRT\n");	
		return 0;
	}


	printf("End reading chemical inputfile\n");
	return 1;
}

int solver::readParam(){

	FILE *paramfile;
	paramfile = fopen("param.dat", "r");

	if(paramfile == NULL){
		printf("Error, param.dat file does not exist\n");
		return 0;
	}


	//------------------------------------------------
	//set default values
	sprintf(inputFilename, "-");
	sprintf(outputFilename, "output.dat");
	//------------------------------------------------

	char sp[160];
	int er; 

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the param.dat file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(paramfile);
			if(c == EOF){
				break;
			}
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
		}
		if(c == EOF) break;
		if(strcmp(sp, "Nwalker =") == 0){
			er = fscanf (paramfile, "%d", &N);
			if(er <= 0){
				printf("Error: Nwalker is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "nSteps =") == 0){
			er = fscanf (paramfile, "%d", &nSteps);
			if(er <= 0){
				printf("Error: nSteps is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "nIterations =") == 0){
			er = fscanf (paramfile, "%d", &nIterations);
			if(er <= 0){
				printf("Error: nIterations is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "nBestSolutions =") == 0){
			er = fscanf (paramfile, "%d", &nBestSolutions);
			if(er <= 0){
				printf("Error: nBestSolutions is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "outputInterval =") == 0){
			er = fscanf (paramfile, "%d", &outputInterval);
			if(er <= 0){
				printf("Error: outputInterval is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "method =") == 0){
			er = fscanf (paramfile, "%d", &method);
			if(er <= 0){
				printf("Error: method is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "Use svgd =") == 0){
			er = fscanf (paramfile, "%d", &svgd);
			if(er <= 0){
				printf("Error: Use svgd is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "f1 =") == 0){
			er = fscanf (paramfile, "%lf", &f1);
			if(er <= 0){
				printf("Error: f1 is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "f2 =") == 0){
			er = fscanf (paramfile, "%lf", &f2);
			if(er <= 0){
				printf("Error: f2 is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "eta =") == 0){
			er = fscanf (paramfile, "%lf", &initial_eta);
			if(er <= 0){
				printf("Error: eta is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "lambda_reg =") == 0){
			er = fscanf (paramfile, "%lf", &lambda_reg);
			if(er <= 0){
				printf("Error: lambda_reg is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "Initial conditions file =") == 0){
			er = fscanf (paramfile, "%s", inputFilename);
			if(er <= 0){
				printf("Error: Initial conditions file is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		if(strcmp(sp, "Output file =") == 0){
			er = fscanf (paramfile, "%s", outputFilename);
			if(er <= 0){
				printf("Error: Output file is not valid!\n");
				return 0;
			}
			if(fgets(sp, 3, paramfile) != nullptr)
			continue;
		}
		printf("Error: param.dat file is not valid! %s\n", sp);
		return 0;
	}
	fclose(paramfile);

	if(strcmp(inputFilename, "-") != 0){
		useICFile = 1;
	} 

	return 1;
}

//read IC file for a single walker
int solver::readIC(){

	printf("Start reading initial conditions file %s\n", inputFilename);

	FILE *file;
	file = fopen(inputFilename, "r");

	if(file == NULL){
		printf("Error, %s file does not exist\n", inputFilename);
		return 0;
	}
 
	int er = 0;
	char sp[160];
	char st[160];
	int pCount = 0;

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the initial conditions file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(file);
			if(c == EOF){
				break;
			}
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
			if(c == '\n'){
				//blank line
				i = -1;
				continue;
			}
		}
		if(c == EOF) break;
//printf("j = %d %s\n", j, sp);


		//Read initial values array
		int fp = 0;
		for(int p = 0; p < Ndim; ++p){
			sprintf(st, "%s =", VariableNames[p].c_str());
			if(strcmp(sp, st) == 0){
				er = fscanf (file, "%lf", &initial_values[p]);
//printf("%s %g | %d\n", st, initial_values[p], p);
				if(er <= 0){
					printf("Error: %s is not valid!\n", VariableNames[p].c_str());
					return 0;
				}
				if(fgets(sp, 3, file) != nullptr)
				fp = 1;
				++pCount;
			}
		}
		if(fp == 1) continue;



		printf("Error in reading initial.dat file\n");
		return 0;
	}

	fclose(file);

	if(pCount != Ndim){
		printf("Error, number of variables in %s file does not agree with Ndim\n", inputFilename);
		return 0;
	}


	printf("End reading initial conditions file %s\n", inputFilename);
	return 1;
}

