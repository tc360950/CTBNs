#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H
#include <vector>

#include "../chain_structure/transition_repository.h"
#include "../chain_structure/state.h"

template <class Real_t> class LikelihoodCalculator {
private:
	TransitionRepository<Real_t> transition_repository;

	Real_t get_intensity(const std::vector<Real_t> &beta, const std::vector<Real_t> &predictive) const {
		Real_t result = 0;
		for (size_t i = 0; i < predictive.size(); i++) {
				result += predictive[i] * beta[i];
		}
		return result;
	}

public:
	LikelihoodCalculator<Real_t>(TransitionRepository<Real_t> transition_repository, bool add_interactions) : 
		transition_repository{ transition_repository } {
	}

	size_t get_parameters_size() const {
		return transition_repository.get_parameters_size();
	}

	Real_t calculate_likelihood(const std::vector<Real_t> &beta, size_t node, bool past_node_value) {
		Real_t result = 0.0;
		const NodeTransitions<Real_t> &node_transitions = transition_repository.fetch_node_transitions(node, past_node_value);
		for (size_t i = 0; i < node_transitions.state_counts.size(); i++) {
			const Real_t intensity = get_intensity(beta, node_transitions.predictive_vectors[i]);
			const Real_t multiplier = -node_transitions.state_counts[i] * intensity 
				+ std::exp(intensity) * node_transitions.time_spent_in_state[i];
			result += multiplier;
		}
		return result;
	}

	std::vector<Real_t> calculate_likelihood_gradient(const std::vector<Real_t> &beta, size_t node, size_t past_node_value) {
		std::vector<Real_t> result(beta.size());
		const NodeTransitions<Real_t> &node_transitions = transition_repository.fetch_node_transitions(node, past_node_value);
		for (size_t i = 0; i < node_transitions.state_counts.size(); i++) {
			const Real_t intensity = get_intensity(beta, node_transitions.predictive_vectors[i]);
			const Real_t multiplier = -node_transitions.state_counts[i] + std::exp(intensity) * node_transitions.time_spent_in_state[i];
			for (size_t j = 0; j < node_transitions.predictive_vectors[i].size(); j++) {
				result[j] += node_transitions.predictive_vectors[i][j] * multiplier;
			}
		}
		return result;
	}
};


#endif // !LIKELIHOOD_CALCULATOR_H
