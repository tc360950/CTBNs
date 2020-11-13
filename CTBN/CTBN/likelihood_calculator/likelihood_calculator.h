#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H

#include "../chain_structure/transition_repository.h"

template <class Real_t> class LikelihoodCalculator {
private:
	TransitionRepository<Real_t> transition_repository;
	const boolean ADD_INTERACTIONS;

	void convert_to_predictive_vector(const State &state, size_t node, std::vector<bool> &result_placeholder) {
		if (!ADD_INTERACTIONS) {
			for (size_t i = 0; i < node; i++) {
				result_placeholder[i] = state.get_node_value(i);
			}
			for (size_t i = node; i < result_placeholder.size(); i++) {
				result_placeholder[i] = state.get_node_value(i + 1);
			}
		}
		return nullptr;
	}

	Real_t get_intensity(const std::vector<Real_t> &beta, const std::vector<bool> &predictive) {
		Real_t result = 0;
		for (size_t i = 0; i < predictive.size(); i++) {
			if (predictive[i) {
				result += beta[i];
			}
		}
		return result;
	}

public:
	LikelihoodCalculator<Real_t>(TransitionRepository<Real_t> transition_repository, boolean add_interactions) : 
		transition_repository{ transition_repository },
		ADD_INTERACTIONS{ add_interactions } {

	}

	Real_t calculate_likelihood(const std::vector<Real_t> &beta, size_t node, bool past_node_value) const {
		Real_t result = 0.0;
		std::vector<bool> predictive_vector(beta.size());
		for (auto &transition_count : transition_repsitory.fetch_node_transitions(node).get_transition_counts()) {
			Transition &transition = transition_count.first;
			if (transition.get_state()[node] == past_node_value) {
				convert_to_predictive_vector(transition.get_state(), node, predictive_vector);
				const Real_t intensity = get_intensity(beta, predictive_vector);
				result += -transition_count.second * intensity + std::exp(intensity) * transition_repository.get_occupation_time(transition.get_state());
			}
		}
		return result;
	}

	std::vector<Real_t> calculate_likelihood_gradient(const std::vector<Real_t> &beta, size_t node, bool past_node_value) const {
		std::vector<Real_t> result(beta.size());
		std::vector<bool> predictive_vector(beta.size());
		for (auto &transition_count : transition_repsitory.fetch_node_transitions(node).get_transition_counts()) {
			Transition &transition = transition_count.first;
			if (transition.get_state()[node] == past_node_value) {
				convert_to_predictive_vector(transition.get_state(), node, predictive_vector);
				const Real_t intensity = get_intensity(beta, predictive_vector, node));
				const Real_t multiplier = -transition_count.second + std::exp(intensity) * transition_repository.get_occupation_time(transition.get_state());
				for (size_t i = 0; i < predictive_vector.size(); i++) {
					if (predictive_vector[i]) {
						result[i] += multiplier;
					}
				}
			}
		}
		return result;
	}
};


#endif // !LIKELIHOOD_CALCULATOR_H
