#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H

#include "../likelihood_calculator/likelihood_calculator.h"
template <class Real_t> class ADMMSolver {
private:
	LikelihoodCalculator<Real_t> likelihood_calculator;
	const Real_t lambda;
	const Real_t RO = 1.0;
	const size_t L2_ITERATIONS = 100;
	const Real_t L2_STEP_SIZE = 0.1;

	void update_x(const std::vector<Real_t> &u, const std::vector<Real_t> &z, std::vector<Real_t> &x, const size_t node, const bool past_node_value) const {
		for (size_t i = 0; i < L2_ITERATIONS; i++) {
			auto gradient_1 = likelihood_calculator.calculate_likelihood_gradient(x, node, past_node_value);
			for (size_t n = 0; n < gradient_1.size(); n++) {
				gradient_1[n] += RO * (x[n] - (z[n] - u[n]));
				x[n] -= L2_STEP_SIZE * gradient_1[n];
			}
		}
	}

	void update_z(const std::vector<Real_t> &u, const std::vector<Real_t> &x, std::vector<Real_t> &z) const {
		const Real_t k = lambda / RO;
		for (size_t i = 0; i < z.size(); i++) {
			if (x[i] + u[i] > k) {
				z[i] = x[i] + u[i] - k;
			}
			else if (x[i] + u[i] < -k) {
				z[i] = x[i] + u[i] + k;
			}
			else {
				z[i] = 0.0;
			}
		}
	}

	void update_u(const std::vector<Real_t> &z, const std::vector<Real_t> &x, std::vector<Real_t> &u) const {
		for (size_t i = 0; i < u.size(); i++) {
			u[i] = u[i] + x[i] - z[i];
		}
	}
public:
	//TODO incijalizacja wektorow u z
	std::vector<Real_t> solve(const size_t iterations, const std::vector<Real_t> &starting_beta, const size_t node, const bool past_node_value) const {
		std::vector<Real_t> x(starting_beta.begin(), starting_beta.end());
		std::vector<Real_t> z;
		std::vector<Real_t> u;
		for (size_t i = 0; i < iterations; i++) {
			update_x(u, z, x, node, past_node_value);
			update_z(u, x, z);
			update_u(z, x, u);
		}
		return x;
	}
};
#endif // !ADMM_SOLVER_H
