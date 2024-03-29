#ifndef S_COOR_H
#define S_COOR_H
#include <thread>
#include <tuple>
#include <string>

#include "solvers/admm_solver.h"
#include "models/model_data.h"
#include "utils/result.h"
#include "utils/parameters.h"

template <class Real_t, class Model> class SamplingCoordinator {
private:
	const size_t SOLVER_ITERATIONS = 10;
	Real_t time_max_local;

	Real_t prune(const std::vector<Real_t> &vector, const Real_t delta, std::vector<Real_t> &result_place_holder) const {
		Real_t non_zero_entries = 0.0;
		for (size_t i = 0; i < vector.size(); i++) {
			if (std::abs(vector[i]) > delta) {
				result_place_holder[i] = vector[i];
				non_zero_entries++;
			}
			else {
				result_place_holder[i] = 0.0;
			}
		}
		return non_zero_entries;
	}

	std::vector<Real_t> get_lambda_sequence(const Real_t max_lambda) const {
		Real_t log_max_lambda = std::log(max_lambda);
		Real_t log_min_lambda = std::log(max_lambda) + std::log(0.0001);
		Real_t diff = (log_max_lambda - log_min_lambda) / 100;
		std::vector<Real_t> result; 
		result.push_back(max_lambda);
		Real_t curr_lambda = log_max_lambda;
		for (size_t i = 0; i < 99; i++) {
			curr_lambda -= diff;
			result.push_back(std::exp(curr_lambda));
		}
		return result;
	}


	Result<Real_t> solve(const Real_t number_of_nodes, ADMMSolver<Real_t> &solver, const Real_t total_number_of_jumps,
		const size_t node, const size_t past_node_value, const TransitionRepository<Real_t> transition_repository) const {
		Real_t number_of_jumps = total_number_of_jumps;
		if (N_DEFINITION == 1) {
			number_of_jumps = transition_repository.fetch_node_transitions(node, past_node_value).number_of_jumps_from_node_value;
		}
		else if (N_DEFINITION == 2) {
			number_of_jumps = transition_repository.fetch_node_transitions(node, past_node_value).number_of_value_changing_jumps;
		}
		std::vector<Real_t> lambdas = get_lambda_sequence(solver.get_max_lambda(node, past_node_value));
	
		std::vector<Real_t> best_so_far;
		Real_t best_so_far_score = 0.0;
		bool best_set = false;
		for (size_t i = 0; i < lambdas.size(); i++) {
			auto result = solver.solve(SOLVER_ITERATIONS, node, past_node_value, lambdas[i]);
			Real_t score = std::get<1>(result) + std::log(number_of_jumps) * std::get<2>(result) / time_max_local;
			if (N_DEFINITION == 3) {
				score = time_max_local * std::get<1>(result) + std::log(time_max_local) * std::get<2>(result);
			}
			if (!best_set || score < best_so_far_score) {
				best_set = true;
				best_so_far_score = score;
				best_so_far = std::get<0>(result);
			}
		}
		if (DEBUG) {
            log("Choosing best delta for node, value pair: ", node, " ", past_node_value);
            log("Chosen vector: "); 
            std::string vector;
            for (auto v : best_so_far) {
                vector.append(std::to_string(v));
                vector.append(" ");
            }
            log(vector);
        } 
		std::vector<Real_t> delta_sequence;
		for (auto el : best_so_far) {
			delta_sequence.push_back(std::abs(el));
		}
		std::vector<Real_t> prunning_place_holder = best_so_far;
		best_set = false;
		std::vector<Real_t> best_prunned_so_far;
        Real_t best_delta;
		for (auto delta : delta_sequence) {
			Real_t non_zero_entries = prune(best_so_far, delta, prunning_place_holder);
			Real_t score =  solver.likelihood_calculator.calculate_likelihood(prunning_place_holder, node, past_node_value);
			score += std::log(2 * number_of_nodes * (number_of_nodes - 1)) * non_zero_entries / time_max_local;
			if (N_DEFINITION == 3) {
				score = time_max_local * solver.likelihood_calculator.calculate_likelihood(prunning_place_holder, node, past_node_value) 
							+ std::log(2 * number_of_nodes * (number_of_nodes - 1)) * non_zero_entries;
			}
            if (DEBUG) {
                log("Delta : ", delta, " total score: ", score, " number fo jumps: ", number_of_jumps);
                log("Likelihood: ",  solver.likelihood_calculator.calculate_likelihood(prunning_place_holder, node, past_node_value), " non_zero_entries ", non_zero_entries);
            }
			if (!best_set || score < best_so_far_score) {
				best_set = true;
				best_so_far_score = score;
				best_prunned_so_far = prunning_place_holder;
                best_delta = delta;
			}
			prunning_place_holder = best_so_far;
		}
        if (DEBUG) {
            log("Best delta: ", best_delta);
        }
		return Result<Real_t>{ best_prunned_so_far, node, past_node_value };
	}

	std::pair<ADMMSolver<Real_t>, ModelData<Real_t>> sample_chain(const size_t node_count, const long seed, const Real_t t_max) const {
		Model model{ node_count, seed };
		auto model_data = model.sample_chain(t_max);
		LikelihoodCalculator<Real_t> calculator{model_data.transition_repository, t_max};
		return std::make_pair(ADMMSolver<Real_t>(calculator), model_data);
	}
public:
	SamplingCoordinator<Real_t, Model>(){}

	SimulationResult<Real_t> simulate(const size_t node_count, const long seed, const Real_t t_max) {
		std::vector<Result<Real_t>> inference_result(2 * node_count);
		auto chain = sample_chain(node_count, seed, t_max);
		time_max_local = t_max;
		std::vector<std::thread> threads;
		const size_t data_per_thread = 2 * node_count / NUM_THREADS;
		for (int th = 0; th < NUM_THREADS; th++) {
			threads.emplace_back([this, th, data_per_thread, node_count, &chain, &inference_result] { 
				const size_t start = data_per_thread * th;
				const size_t end = th + 1 == NUM_THREADS ? 2 * node_count : start + data_per_thread;
				for (size_t i = start; i < end; i++) {
					const size_t node_value = i % 2;
					auto result = this->solve(node_count, chain.first, chain.second.number_of_jumps, i / 2, node_value, chain.second.transition_repository);
					inference_result[i] = result;
				}
			});
		}
		for (auto &th : threads)
		{
			th.join();
		}
		return SimulationResult<Real_t>{ inference_result, chain.second };
	}

};
#endif // !S_COOR_H_H
