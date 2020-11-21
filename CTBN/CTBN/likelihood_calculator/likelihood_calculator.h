#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H
#include <vector>
#include <algorithm>

#include "../chain_structure/transition_repository.h"
#include "../chain_structure/state.h"
#include "../parameters/parameters.h"

template <class Real_t> class LikelihoodCalculator {
private:
	TransitionRepository<Real_t> transition_repository;
	const Real_t t_max;

	inline Real_t get_intensity(const std::vector<Real_t> &beta, const std::vector<Real_t> &predictive, const size_t start) const {
		Real_t result = 0;
		for (size_t i = 0; i < beta.size(); i++) {
			result += predictive[start + i] * beta[i];
		}
		return result;
	}

public:
	LikelihoodCalculator<Real_t>(TransitionRepository<Real_t> transition_repository, Real_t t_max) :
		transition_repository{ transition_repository },
		t_max{ t_max } {
	}

	size_t get_parameters_size() const {
		return transition_repository.get_parameters_size();
	}

	Real_t calculate_likelihood(const std::vector<Real_t> &beta, size_t node, bool past_node_value) {
		Real_t result = 0.0;
		size_t start_od_predictive = 0;
		const NodeTransitions<Real_t> &node_transitions = transition_repository.fetch_node_transitions(node, past_node_value);
		for (size_t i = 0; i < node_transitions.state_counts.size(); i++) {
			const Real_t intensity = get_intensity(beta, node_transitions.predictive_vectors, start_od_predictive);
			const Real_t multiplier = -node_transitions.state_counts[i] * intensity
				+ std::exp(intensity) * node_transitions.time_spent_in_state[i];
			result += multiplier;
			start_od_predictive += get_parameters_size();
		}
		return result / t_max;
	}

	void calculate_likelihood_gradient(const std::vector<Real_t> &beta, size_t node, size_t past_node_value, std::vector<Real_t> &result) {
		size_t start_od_predictive = 0;
		std::fill(result.begin(), result.end(), 0.0);
		const NodeTransitions<Real_t> &node_transitions = transition_repository.fetch_node_transitions(node, past_node_value);
		for (size_t i = 0; i < node_transitions.state_counts.size(); i++) {
			const Real_t exp_intensity = std::exp(get_intensity(beta, node_transitions.predictive_vectors, start_od_predictive));
			for (size_t j = 0; j < beta.size(); j++) {
				result[j] += node_transitions.predictive_vectors_times_occupation[start_od_predictive + j] * exp_intensity;
			}
			start_od_predictive += get_parameters_size();
		}
		std::transform(result.begin(), result.end(), node_transitions.sum_counts_times_predictive.begin(),
			result.begin(), std::plus<Real_t>());
		for (size_t j = 0; j < beta.size(); j++) {
			result[j] = result[j] / t_max;
		}
	}
};


#endif // !LIKELIHOOD_CALCULATOR_H