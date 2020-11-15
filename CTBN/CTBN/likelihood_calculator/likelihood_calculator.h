#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H
#include <vector>

#include "../chain_structure/transition_repository.h"
#include "../chain_structure/state.h"

template <class Real_t> class LikelihoodCalculator {
private:
	TransitionRepository<Real_t> transition_repository;
	std::vector<std::unordered_map<State, std::vector<Real_t>, StateHash>> predictive_vectors_cache;
	const bool ADD_INTERACTIONS;

	const std::vector<Real_t> &convert_to_predictive_vector(const State &state, size_t node) {
		if (!ADD_INTERACTIONS) {
			if (predictive_vectors_cache[node].find(state) != predictive_vectors_cache[node].end()) {
				return predictive_vectors_cache[node][state];
			}
			std::vector<Real_t> result;
			for (size_t i = 0; i < node; i++) {
				result.push_back(state.get_node_value(i));
			}
			for (size_t i = node; i < state.get_size() - 1; i++) {
				result.push_back(state.get_node_value(i + 1));
			}
			predictive_vectors_cache[node][state] = result;
			return predictive_vectors_cache[node][state];
		}
	}

	Real_t get_intensity(const std::vector<Real_t> &beta, const std::vector<Real_t> &predictive) const {
		Real_t result = 0;
		for (size_t i = 0; i < predictive.size(); i++) {
				result += predictive[i] * beta[i];
		}
		return result;
	}

public:
	LikelihoodCalculator<Real_t>(TransitionRepository<Real_t> transition_repository, bool add_interactions) : 
		transition_repository{ transition_repository },
		ADD_INTERACTIONS{ add_interactions } {
		predictive_vectors_cache.resize(transition_repository.get_number_of_nodes());
	}

	size_t get_parameters_size() const {
		if (!ADD_INTERACTIONS) {
			return transition_repository.get_number_of_nodes() - 1;
		}
	}

	Real_t calculate_likelihood(const std::vector<Real_t> &beta, size_t node, bool past_node_value) {
		Real_t result = 0.0;
		std::vector<Real_t> *predictive_vector;
		for (auto &transition_count : transition_repository.fetch_node_transitions(node).get_transition_counts()) {
			const State &state = transition_count.first.get_state();
			if (state.get_node_value(node) == past_node_value) {
				const std::vector<Real_t> &predictive_vector = convert_to_predictive_vector(state, node);
				const Real_t intensity = get_intensity(beta, predictive_vector);
				result += -transition_count.second * intensity + std::exp(intensity) * transition_repository.get_occupation_time(state);
			}
		}
		return result;
	}

	std::vector<Real_t> calculate_likelihood_gradient(const std::vector<Real_t> &beta, size_t node, bool past_node_value) {
		std::vector<Real_t> result(beta.size());
		for (auto &transition_count : transition_repository.fetch_node_transitions(node).get_transition_counts()) {
			const State &state = transition_count.first.get_state();
			if (state.get_node_value(node) == past_node_value) {
				const std::vector<Real_t> &predictive_vector = convert_to_predictive_vector(state, node);
				const Real_t intensity = get_intensity(beta, predictive_vector);
				const Real_t multiplier = -transition_count.second + std::exp(intensity) * transition_repository.get_occupation_time(state);
				for (size_t i = 0; i < predictive_vector.size(); i++) {
						result[i] += predictive_vector[i] * multiplier;
				}
			}
		}
		return result;
	}
};


#endif // !LIKELIHOOD_CALCULATOR_H
