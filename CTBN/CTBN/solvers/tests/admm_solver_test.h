#ifndef ADMM_SOLVER_TEST_H
#define ADMM_SOLVER_TEST_H
#include <vector>

#include "../admm_solver.h"
#include "../../utils/logger.h"

template <class Real_t, class Model> class ADMMSolverTest {

private:
	const size_t NUMBER_OF_NODES = 20;
	const long SEED = 213141;
	const Real_t T_MAX = 50.0;
	const Real_t TOLERANCE = 0.01;
	const size_t NUMBER_OF_TESTS = 5;
	std::mt19937 generator{ SEED };

	std::vector<Real_t> get_starting_x(const size_t vector_size) {
		std::uniform_real_distribution<Real_t> distribution(0.0, 1.0);
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(distribution(generator));
		}
		return result;
	}

	std::vector<Real_t> get_starting_z(const size_t vector_size) {
		std::uniform_real_distribution<Real_t> distribution(0.0, 1.0);
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(distribution(generator));
		}
		return result;
	}

	std::vector<Real_t> get_starting_u(const size_t vector_size) {
		std::uniform_real_distribution<Real_t> distribution(0.0, 1.0);
		std::vector<Real_t> result;
		for (size_t i = 0; i < vector_size; i++) {
			result.push_back(distribution(generator));
		}
		return result;
	}

	std::vector<Real_t> solve(LikelihoodCalculator<Real_t> calculator, const size_t node, const size_t past_node_value, const Real_t lambda) {
		ADMMSolver<Real_t> admm_solver{ calculator };
		std::vector<Real_t> x = get_starting_x(NUMBER_OF_NODES);
		std::vector<Real_t> z = get_starting_z(NUMBER_OF_NODES);
		std::vector<Real_t> u = get_starting_u(NUMBER_OF_NODES);
		std::vector<Real_t> gradient_holder;
		gradient_holder.resize(x.size());
		bool stop = false;
		size_t counter = 0;
		while (!stop) {
			for (size_t i = 0; i < 10 - 1; i++) {
				admm_solver.update_x(u, z, x, node, past_node_value, gradient_holder);
				admm_solver.update_z(u, x, z, lambda);
				admm_solver.update_u(z, x, u);
				counter++;
			}
			std::vector<Real_t> previous_u = u;
			std::vector<Real_t> previous_z = z;
			admm_solver.update_x(u, z, x, node, past_node_value, gradient_holder);
			admm_solver.update_z(u, x, z, lambda);
			admm_solver.update_u(z, x, u);
			stop = admm_solver.check_stop(u, x, z, previous_u, previous_z);
			counter++;
		}
		return z;
	}

public:
	void test(long seed, Real_t lambda) {
		Model model{ NUMBER_OF_NODES, seed };
		auto model_data = model.sample_chain_and_skeleton(T_MAX);
		LikelihoodCalculator<Real_t> calculator{ model_data.first.transition_repository, T_MAX };
		for (size_t i = 0; i < 2 * NUMBER_OF_NODES; i++) {
			auto node = i / 2;
			auto value = i % 2;
			auto z = solve(calculator, node, value, lambda);
			bool ok = true;
			for (size_t t = 0; t < NUMBER_OF_TESTS; t++) {
				std::cout << "ADMM solver test start\n";
				auto z_2 = solve(calculator, node, value, lambda);
				std::cout << "ADMM solver test end\n";
				for (size_t i = 0; i < z.size(); i++) {
					if (std::abs(z[i] - z_2[i]) >= TOLERANCE) {
						logTest<ADMMSolverTest>("Test fail at ", z[i], "vs", z_2[i]);
						while (true) { std::cout << "dupa"; };
						ok = false;
					}
				}
			}
			if (ok) {
				logTest<ADMMSolverTest>("Test ok!");
			}
		}
	}
};
#endif // !ADMM_SOLVER_TEST_H
