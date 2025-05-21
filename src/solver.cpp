#include "solver.h"



//Allocate all arrays
//N is the number of chains,
//Ndim is the number of equations
void solver::allocate(){

	GRT_T = (double*)malloc((GRTmaxIndex + 1) * sizeof(double));
	prior_min = (double*)malloc(Ndim * sizeof(double));
	prior_max = (double*)malloc(Ndim * sizeof(double));
	bound_min = (double*)malloc(Ndim * sizeof(double));
	bound_max = (double*)malloc(Ndim * sizeof(double));
	logbound_min = (double*)malloc(Ndim * sizeof(double));
	logbound_max = (double*)malloc(Ndim * sizeof(double));

	//set prior and boundaries to a negative value to check
	//if they are given in the chem_input file.
	for(int i = 0; i < Ndim; ++i){
		prior_min[i] = LARGE;
		prior_max[i] = LARGE;
		bound_min[i] = LARGE;
		bound_max[i] = LARGE;
	}
	for(int i = 0; i < GRTmaxIndex + 1; ++i){
		GRT_T[i] = LARGE;
	}

	for(int i = 0; i < nParameters; ++i){
		Parameters[i] = LARGE;
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


int solver::readGibbs(){

	if(useGibbsFile == 0){
		printf("Use  Gibbs.py script to generate Gibbs energies file Gibbs.dat\n");
		char command[160];
		sprintf(command, "python3 ../Gibbs.py %g %g", T_surf, T_CMB);
		int er = system(command);
		sprintf(GibbsFilename, "Gibbs.dat");

		if(er != 0){
			printf("Error in running Gibbs.py: %s\n", command);
			return 0;

		}

	}




	printf("Start reading Gibbs energy file %s\n", GibbsFilename);

	FILE *Gibbsfile;
	Gibbsfile = fopen(GibbsFilename, "r");

	if(Gibbsfile == NULL){
		printf("Error, Gibbs file does not exist\n");
		return 0;
	}

	//------------------------------------------------

	char sp[160];
	char st[160];

	int er; 
	int grtCount = 0;

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the Gibbs file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(Gibbsfile);
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
		double T;
		for(int p = 0; p < nGRT; ++p){
			sprintf(st, "GRT_%d =", GRTid[p]);

			if(strcmp(sp, st) == 0){
				int id = GRTid[p];
				er = fscanf (Gibbsfile, "%lf %lf", &GRT_T[id], &T);
//printf("%s %d %d %g %g %g %g\n", sp, p, id, GRT_T[id], T, T_surf, T_CMB);
				if(er <= 0){
					printf("Error: GRT_%d is not valid!\n", p);
					return 0;
				}
				//check if Temperatures agree
				if(abs(T - T_surf) > 5.0 && abs(T - T_CMB) > 5.0){
					printf("Error: Temperature for Gibbs Energy GRT_ %d does not math T_surf or T_CMB. %g %g %g\n", p, T, T_surf, T_CMB);
					return 0;

				}


				if(fgets(sp, 3, Gibbsfile) != nullptr)
				fg = 1;
				++grtCount;
			}
		}
		if(fg == 1) continue;

	}

	fclose(Gibbsfile);

	//Check if all GRT values are set
	printf("T_surf %g T_CMB %g\n", T_surf, T_CMB);
	for(int i = 0; i < nGRT; ++i){
		int id = GRTid[i];
		printf("GRT %d %g\n", id, GRT_T[id]);
		if(GRT_T[id] == LARGE){
			printf("Error, value for GRT %d not set in Gibbs file\n", id);
			return 0;

		}
	}

	if(grtCount != nGRT){
		printf("Error, number of GRT terms in Gibbs file %d does not agree with nGRT %d\n", grtCount, nGRT);	
		return 0;
	}
	printf("End reading Gibbs energy file\n");

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
	int parameterCount = 0;
	int boundariesCount = 0;
	int priorsCount = 0;

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the chem_input.dat file
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
				++parameterCount;
			}
		}
		if(fp == 1) continue;

		//Read Boundaries array
		int fbound = 0;
		for(int p = 0; p < Ndim; ++p){
			sprintf(st, "bound_%s =", VariableNames[p].c_str());
			if(strcmp(sp, st) == 0){
				er = fscanf (paramfile, "%lf %s %lf", &bound_min[p], skip, &bound_max[p]);
//printf("%s %g %g| %d\n", st, bound_min[p], bound_max[p], p);
				logbound_min[p] = log(bound_min[p]);
				logbound_max[p] = log(bound_max[p]);
				if(er <= 0 || bound_min[p] < 0.0 || bound_max[p] < 0.0){
					printf("Error: bound_%s is not valid!\n", VariableNames[p].c_str());
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fbound = 1;
				++boundariesCount;
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
				if(er <= 0 || prior_min[p] < 0.0 || prior_max[p] < 0.0){
					printf("Error: prior_%s is not valid!\n", VariableNames[p].c_str());
					return 0;
				}
				if(fgets(sp, 3, paramfile) != nullptr)
				fprior = 1;
				++priorsCount;
			}
		}
		if(fprior == 1) continue;

		
		printf("Error: chem_input.dat file is not valid! %s\n", sp);
		return 0;
	}
	fclose(paramfile);


	if(parameterCount != nParameters){
		printf("Error, number of parameters in chem_input.dat file does not agree with nParameters\n");	
		return 0;
	}
	if(boundariesCount != Ndim){
		printf("Error, number of bondaries in chem_input.dat file does not agree with number of variables\n");	
		return 0;
	}
	if(priorsCount != Ndim){
		printf("Error, number of priors in chem_input.dat file does not agree with number of variables\n");	
		return 0;
	}


	//Check if all Parameters values are set
	for(int i = 0; i < nParameters; ++i){
		if(Parameters[i] == LARGE){
			printf("Error, value for Parameter %s not set in chem_input.dat file\n", ParameterNames[i].c_str());
			return 0;

		}
	}

	//Check if all priors and boundaries are set
	for(int i = 0; i < Ndim; ++i){
		if(bound_min[i] == LARGE || bound_max[i] == LARGE){

			printf("Error, boudaries for %s not set in chem_input.dat file\n", VariableNames[i].c_str());
			return 0;

		}
		if(prior_min[i] == LARGE || prior_max[i] == LARGE){

			printf("Error, priors for %s not set in chem_input.dat file\n", VariableNames[i].c_str());
			return 0;

		}
	}


	printf("End reading chemical inputfile\n");



	FILE *outfile;
	outfile = fopen(outputFilename, "r"); // Check if output file exists already
	if(outfile == NULL){
		//File does not exist, add header information
		outfile = fopen(outputFilename, "w");
		fprintf(outfile, "#iteration chain chi^2 ");
		for(int i = 0; i < Ndim; ++i){
			fprintf(outfile, "%s ", VariableNames[i].c_str()); 
		}
		fprintf(outfile, "\n");
		fclose(outfile);
	}
	else{
		//File exists already
		char header[160];
		int er;
		for(int i = 0; i < 3; ++i){
			er = fscanf(outfile, "%s", header);
			if(er <= 0){
				printf("Error in reading the header of the output file\n");
				return 0;
			}
		}
		for(int i = 0; i < Ndim; ++i){
			er = fscanf(outfile, "%s", header);
			if(er <= 0){
				printf("Error in reading the header of the output file\n");
				return 0;
			}
			int check = strcmp(header, VariableNames[i].c_str());
			if(check != 0){
				printf("Error, the header of the output file does not match the current variables\n");
				printf("%s %s\n", header, VariableNames[i].c_str());
				return 0;
			}
		}
		fclose(outfile);
	}


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
	sprintf(GibbsFilename, "Gibbs.dat");
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
			if(c == '\n'){
				//blank line
				i = -1;
				continue;
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
		if(strcmp(sp, "Stop condition =") == 0){
			er = fscanf (paramfile, "%lf", &stopcritera);
			if(er <= 0){
				printf("Error: Stop condition is not valid!\n");
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
		if(strcmp(sp, "Gibbs energy file =") == 0){
			er = fscanf (paramfile, "%s", GibbsFilename);
			if(er <= 0){
				printf("Error: Gibbs file is not valid!\n");
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
	if(strcmp(GibbsFilename, "-") != 0){
		useGibbsFile = 1;
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

