#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H
#include <vector>
#include <tuple>

#include "../likelihood_calculator/likelihood_calculator.h"
#include "../bob_dylan.h"
#include "../utils/logger.h"
#include "tests/admm_solver_test.h"

template <class Real_t> class ADMMSolver {
public:
	LikelihoodCalculator<Real_t> likelihood_calculator;
private:
	const Real_t RO = 1.0;
	const size_t L2_ITERATIONS = 10;
	const Real_t L2_STEP_SIZE = 0.1;
	const Real_t EPSILON = 0.000000001;
	const Real_t STOPPING_EPSILON = 0.0001;

	void update_x(const std::vector<Real_t> &u, const std::vector<Real_t> &z, std::vector<Real_t> &x, const size_t node, const bool past_node_value, std::vector<Real_t> &gradient_holder) {
		bool stop = false;
		while (!stop) {
			likelihood_calculator.calculate_likelihood_gradient(x, node, past_node_value, gradient_holder);
			for (size_t n = 0; n < gradient_holder.size(); n++) {
				gradient_holder[n] += RO * (x[n] - (z[n] - u[n]));
				x[n] -= L2_STEP_SIZE * gradient_holder[n];
			}
			auto gradient_norm = get_vector_l2_norm(gradient_holder);
			if (gradient_norm <= STOPPING_EPSILON) {
				stop = true;
			}
		}
	}

	void update_z(const std::vector<Real_t> &u, const std::vector<Real_t> &x, std::vector<Real_t> &z, const Real_t lambda) {
		const Real_t k = lambda / RO;
		z[0] = x[0] + u[0];
		for (size_t i = 1; i < z.size(); i++) {
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
		for (size_t i = 1; i < vec.size(); i++) {
			result += std::abs(vec[i]);
		}
		return result;
	}

	Real_t get_vector_l2_norm(const std::vector<Real_t> &vec) const {
		Real_t result = 0.0;
		for (size_t i = 0; i < vec.size(); i++) {
			result += vec[i] * vec[i];
		}
		return std::pow(result, 0.5);
	}

	Real_t difference_norm(const std::vector<Real_t> &vec1, const std::vector<Real_t> &vec2) const {
		Real_t result = 0.0;
		for (size_t i = 0; i < vec1.size(); i++) {
			result += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
		}
		return std::pow(result, 0.5);
	}

	bool check_stop(const std::vector<Real_t> &u, const std::vector<Real_t> &x, const std::vector<Real_t> &z, const std::vector<Real_t> &old_u,
		const std::vector<Real_t> &old_z) const {
		Real_t e_primal = std::sqrt(x.size()) * STOPPING_EPSILON + STOPPING_EPSILON * std::max(get_vector_l2_norm(x), get_vector_l2_norm(z));
		Real_t e_dual = std::sqrt(x.size()) * STOPPING_EPSILON;
		return e_primal >= difference_norm(u, old_u) && e_dual >= RO * difference_norm(z, old_z);
	}
public:
	ADMMSolver<Real_t>(LikelihoodCalculator<Real_t> likelihood_calculator) :
		likelihood_calculator{ likelihood_calculator } {}
	// pierwszy element return to beta, drugi n * lik trzeci ||beta||
	std::tuple<std::vector<Real_t>, Real_t, Real_t> solve(const size_t iterations, const size_t node, const bool past_node_value, const Real_t lambda) {
		std::vector<Real_t> x = get_starting_x(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> z = get_starting_z(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> u = get_starting_u(likelihood_calculator.get_parameters_size());
		std::vector<Real_t> gradient_holder;
		gradient_holder.resize(x.size());
		if (DEBUG) {
			log("Starting ADMM");
			log("Score: ", get_vector_penalty(z) * lambda + likelihood_calculator.calculate_likelihood(z, node, past_node_value));
			log("Penalty: ", get_vector_penalty(z) * lambda, " Likelihood: ", likelihood_calculator.calculate_likelihood(z, node, past_node_value));
		}
		bool stop = false;
		size_t counter = 0;
		while (!stop) {
			for (size_t i = 0; i < iterations - 1; i++) {
				update_x(u, z, x, node, past_node_value, gradient_holder);
				update_z(u, x, z, lambda);
				update_u(z, x, u);
				counter++;
			}
			std::vector<Real_t> previous_u = u;
			std::vector<Real_t> previous_z = z;
			update_x(u, z, x, node, past_node_value, gradient_holder);
			update_z(u, x, z, lambda);
			update_u(z, x, u);
			stop = check_stop(u, x, z, previous_u, previous_z);
			counter++;
		}
		if (DEBUG) {
			log("Stopped ADMM after ", counter, " iterations");
			log("Finished optimization for node ", node, " past node value: ", (size_t)past_node_value, " with lambda ", lambda);
			log("Score: ", get_vector_penalty(z) * lambda + likelihood_calculator.calculate_likelihood(z, node, past_node_value));
			log("Penalty: ", get_vector_penalty(z) * lambda, " Likelihood: ", likelihood_calculator.calculate_likelihood(z, node, past_node_value));
			log("Non zero vector elements: ", get_non_zero_vector_elements(z));
		}
		return std::make_tuple(z, likelihood_calculator.calculate_likelihood(z, node, past_node_value), get_non_zero_vector_elements(z));
	}

	friend class ADMMSolverTest<Real_t>;
};
#endif // !ADMM_SOLVER_H
