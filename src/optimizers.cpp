#include "solver.h"


void solver::setBoundaries(double *vars){
	for (int k = 0; k < N; ++k) {
		for (int i = 0; i < Ndim; ++i) {
			int idx = k * Ndim + i;
			if(vars[idx] < bound_min[i]) {
				vars[idx] = bound_min[i];
			}
			if(vars[idx] > bound_max[i]) {
				vars[idx] = bound_max[i];
			}
		}
	}
}

void solver::setBoundariesLog(double *logvars){
	for (int k = 0; k < N; ++k) {
		for (int i = 0; i < Ndim; ++i) {
			int idx = k * Ndim + i;
			if(logvars[idx] < logbound_min[i]) {
				logvars[idx] = logbound_min[i];
			}
			if(logvars[idx] > logbound_max[i]) {
				logvars[idx] = logbound_max[i];
			}
		}
	}
}

void solver::adamW_update(double *vars) {

	eta *= 0.9999;

	for (int k = 0; k < N; ++k) {
		for (int i = 0; i < Ndim; ++i) {
			int idx = k * Ndim + i;

			m[idx] = beta1 * m[idx] + (1.0 - beta1) * svgd_grad[idx];
			v[idx] = beta2 * v[idx] + (1.0 - beta2) * svgd_grad[idx] * svgd_grad[idx];

			double mm = m[idx] / (1.0 - beta1_t);
			double vv = v[idx] / (1.0 - beta2_t);

			d[idx] = -eta * mm / (sqrt(vv) + eps) - eta * lambda_reg * vars[idx];
			vars[idx] += d[idx];

		}
	}

	beta1_t *= beta1;
	beta2_t *= beta2;
}

void solver::adamaxW_update(double *vars) {

	eta *= 0.9999;

	for (int k = 0; k < N; ++k) {
		//for (int k = N - 1; k < N; ++k) {

		for (int i = 0; i < Ndim; ++i) {
			int idx = k * Ndim + i;

			// Update first moment
			m[idx] = beta1 * m[idx] + (1.0 - beta1) * svgd_grad[idx];

			// Update infinity norm of gradient (maximum absolute gradient)
			v[idx] = fmax(beta2 * v[idx], fabs(svgd_grad[idx]));

			// Compute weight update
			d[idx] = -eta * m[idx] / (v[idx] + eps) - eta * lambda_reg * vars[idx];

			// Apply update to parameter
			vars[idx] += d[idx];

		}
	}
}

void solver::SVGD(double *vars){
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			double dist2 = 0.0;
			for (int d = 0; d < Ndim; ++d) {
				double diff = (vars[i * Ndim + d] - vars[j * Ndim + d]) / h[d];
				dist2 += diff * diff;
			}
			double kij = exp(-dist2);
			for (int d = 0; d < Ndim; ++d) {
				// Gradient contribution
				double diff = (vars[i * Ndim + d] - vars[j * Ndim + d]) / (h[d] * h[d]);
				svgd_grad[i * Ndim + d] += f1 * kij * ctot_grad[j * Ndim + d];

				// Repulsive term
				if (i != j) {
					svgd_grad[i * Ndim + d] -= f2 * 2.0 * kij * diff * 0.25;
				}
			}
		}
	}

	// Normalize SVGD gradient
	for (int k = 0; k < Ndim * N; ++k) {
		svgd_grad[k] /= (double)N;
	}
}

void solver::SVGD_1(double *vars, int i){
	for (int j = 0; j < N; ++j) {
		double dist2 = 0.0;
		for (int d = 0; d < Ndim; ++d) {
			double diff = (vars[i * Ndim + d] - vars[j * Ndim + d]) / h[d];
				dist2 += diff * diff;
		}
		double kij = exp(-dist2);
		for (int d = 0; d < Ndim; ++d) {
			// Gradient contribution
			double diff = (vars[i * Ndim + d] - vars[j * Ndim + d]) / (h[d] * h[d]);
			svgd_grad[i * Ndim + d] += f1 * kij * ctot_grad[j * Ndim + d];

			// Repulsive term
			if (i != j) {
				svgd_grad[i * Ndim + d] -= f2 * 2.0 * kij * diff * 0.25;
			}
		}
	}

	// Normalize SVGD gradient
	//for (int k = 0; k < Ndim * N; ++k) {
	//	svgd_grad[k] /= (double)N;
	//}
}


