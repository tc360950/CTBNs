#ifndef BOB_DYLAN_H
#define BOB_DYLAN_H
#include <thread>

#include "solvers/admm_solver.h"
template <class Real_t, class Model> class BobDylan {
private:
	const size_t NUM_THREADS = 1;
	const size_t SOLVER_ITERATIONS = 10000;
	const Real_t MAX_LAMBDA = 100000;
	const size_t LAMBDA_COUNT = 10;
	//TODO fill it
	const std::vector<Real_t> DELTA_SEQUENCE;

	class Result {
	public:
		std::vector<Real_t> beta;
		size_t node;
		size_t node_past_value;
	};

	Real_t prune(const std::vector<Real_t> &vector, const Real_t delta, std::vector<Real_t> &result_place_holder) {
		Real_t non_zero_entries = 0.0;
		for (size_t i = 0; i < vector.size(); i++) {
			if (vector[i] >= delta || vector[i] <= -delta) {
				result_place_holder[i] = vector[i];
				non_zero_entries++;
			}
			else {
				result_place_holder[i] = 0.0;
			}
		}
		return non_zero_entries;
	}

	Result solve(const Real_t number_of_nodes, const ADMMSolver<Real_t> &solver, const Real_t number_of_jumps, const std::vector<Real_t> &starting_beta, const size_t node, const bool past_node_value) const {
		std::vector<Real_t> best_so_far;
		Real_t best_so_far_score = 0.0;
		bool best_set = false;
		Real_t lambda = MAX_LAMBDA;
		for (size_t i = 0; i < LAMBDA_COUNT; i++) {
			auto result = solver.solve(SOLVER_ITERATIONS, starting_beta, node, past_node_value, lambda);
			Real_t score = number_of_jumps * std::get<1>(result) + std::log(number_of_jumps) * std::get<2>(result);
			if (!best_set || score < best_so_far_score) {
				best_set = true;
				best_so_far_score = score;
				best_so_far = std::get<0>(result);
			}
			lambda = lambda * 0.1;
		}
		//Choose best delta 

		std::vector<Real_t> prunning_place_holder = best_so_far;
		best_set = false;
		std::vector<Real_t> best_prunned_so_far;
		for (auto delta : DELTA_SEQUENCE) {
			Real_t non_zero_entries = prune(best_so_far, delta, prunning_place_holder);
			Real_t score = number_of_jumps * solver.likelihood_calculator.calculate_likelihood(prunning_place_holder, node, past_node_value);
			score += std::log(2 * number_of_nodes / (number_of_nodes - 1)) * non_zero_entries;
			if (!best_set || score < best_so_far_score) {
				best_set = true;
				best_so_far_score = score;
				best_prunned_so_far = prunning_place_holder;
			}
		}
		return Result(best_prunned_so_far, node, past_node_value);
	}

	std::pair<ADMMSolver<Real_t>, Real_t> sample_chain(const size_t node_count, const long seed) {

	}
public:


	void simulate(const size_t node_count, const long seed) {
		std::vector<Result> inference_result(2 * node_count);
		auto chain = sample_chain(node_count, seed);
		std::vector<std::thread> threads;
		const size_t data_per_thread = 2 * node_count / NUM_THREADS;
		for (int th = 0; th < NUM_THREADS; th++) {
			threads.emplace_back([this, th, data_per_thread, node_count, &chain, &inference_result] { 
				const size_t start = data_per_thread * th;
				const size_t end = th + 1 == this->NUM_THREADS ? 2 * node_count : start + data_per_thread;
				for (size_t i = start; i < end; i++) {
					const bool node_value = i % 2 == 1;
					//todo starting beta
					auto result = this->solve(node_count, chain.first, chain.second, starting_beta, i / 2, node_value);
					inference_result[i] = result;
				}
			});
		}
		for (auto &th : threads)
		{
			th.join();
		}
	}

};
#endif // !BOB_DYLAN_H
