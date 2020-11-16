#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H
#include <vector>
#include <tuple>

#include "../likelihood_calculator/likelihood_calculator.h"
#include "../bob_dylan.h"

template <class Real_t> class ADMMSolver {
public:
	LikelihoodCalculator<Real_t> likelihood_calculator;
private:
	const Real_t RO = 1.0;
	const size_t L2_ITERATIONS = 100;
	const Real_t L2_STEP_SIZE = 0.1;
	const Real_t EPSILON = 0.000000001;
	void update_x(const std::vector<Real_t> &u, const std::vector<Real_t> &z, std::vector<Real_t> &x, const size_t node, const bool past_node_value, std::vector<Real_t> &gradient_holder) {
		for (size_t i = 0; i < L2_ITERATIONS; i++) {
			likelihood_calculator.calculate_likelihood_gradient(x, node, past_node_value, gradient_holder);
			for (size_t n = 0; n < gradient_holder.size(); n++) {
				gradient_holder[n] += RO * (x[n] - (z[n] - u[n]));
				x[n] -= L2_STEP_SIZE * gradient_holder[n];
			}
		}
	}

	void update_z(const std::vector<Real_t> &u, const std::vector<Real_t> &x, std::vector<Real_t> &z, const Real_t lambda) {
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

	void update_u(const std::vector<Real_t> &z, const std::vector<Real_t> &x, std::vector<Real_t> &u) {
		for (size_t i = 0; i < u.size(); i++) {
			u[i] = u[i] + x[i] - z[i];
		}
	}

	Real_t get_non_zero_vector_elements(const std::vector<Real_t> &v) {
		Real_t result = 0.0;
		for (auto element : v) {
			if (element > EPSILON || element < -EPSILON) {
				result += 1.0;
			}
		}
		return result;
	}

	std::vector<Real_t> get_starting_x(const size_t vector_size) {
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(.01);
		}
		return result;
	}

	std::vector<Real_t> get_starting_z(const size_t vector_size) {
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(.01);
		}
		return result;
	}

	std::vector<Real_t> get_starting_u(const size_t vector_size) {
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(.01);
		}
		return result;
	}

	Real_t get_vector_penalty(const std::vector<Real_t> &vec) {
		Real_t result = 0.0;
		for (size_t i = 0; i < vec.size(); i++) {
			result += std::abs(vec[i]);
		}
		return result;
	}
public:
	ADMMSolver<Real_t>(LikelihoodCalculator<Real_t> likelihood_calculator) :
		likelihood_calculator{ likelihood_calculator } {}
	//TODO incijalizacja wektorow u z
	// pierwszy element return to beta, drugi n * lik trzeci ||beta||
	std::tuple<std::vector<Real_t>, Real_t, Real_t> solve(const size_t iterations, const size_t node, const bool past_node_value, const Real_t lambda) {
		std::vector<Real_t> x = get_starting_x(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> z = get_starting_z(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> u = get_starting_u(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> gradient_holder;
		gradient_holder.resize(x.size());
		for (size_t i = 0; i < iterations; i++) {
			update_x(u, z, x, node, past_node_value, gradient_holder);
			update_z(u, x, z, lambda);
			update_u(z, x, u);
			//std::cout << "Score: " << get_vector_penalty(x) * lambda + likelihood_calculator.calculate_likelihood(x, node, past_node_value) << "\n";
		}
		return std::make_tuple(x, likelihood_calculator.calculate_likelihood(x, node, past_node_value), get_non_zero_vector_elements(x));
	}
};
#endif // !ADMM_SOLVER_H
