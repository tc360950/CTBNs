#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H

#include "../chain_structure/transition_repository.h"

template <class Real_t> class LikelihoodCalculator {
private:
	TransitionRepository<Real_t> transition_repository;

	//TODO - konwersja stanu do Z - prosta optynalziacja wyrzuc na zewnatrz tej funkcji tworenie tego wektora i zmieniaj tylko jego wartosci, ale 
	// uwazaj na race przy wielu watkach. - to raczej powinno byc lokalne w liczeniu likelihoodu
	std::vector<bool> convert_to_predictive_vector(const State &state, size_t node) {
		return nullptr;
	}
	//ilcozy nskalrany
	Real_t get_intensity(const std::vector<Real_t> &beta, const std::vector<bool> &predictive) {
		return 0.0;
	}

public:
	Real_t calculate_likelihood(const std::vector<Real_t> &beta, size_t node, bool past_node_value) const {
		Real_t result = 0.0;
		for (auto &transition_count : transition_repsitory.fetch_node_transitions(node).get_transition_counts()) {
			Transition &transition = transition_count.first;
			if (transition.get_state()[node] == past_node_value) {
				const Real_t intensity = get_intensity(beta, convert_to_predictive_vector(transition.get_state(), node));
				result += -transition_count.second * intensity + std::exp(intensity) * transition_repository.get_occupation_time(transition.get_state());
			}
		}
		return result;
	}

	std::vector<Real_t> calculate_likelihood_gradient(const std::vector<Real_t> &beta, size_t node, bool past_node_value) const {
		std::vector<Real_t> result(beta.size());
		for (auto &transition_count : transition_repsitory.fetch_node_transitions(node).get_transition_counts()) {
			Transition &transition = transition_count.first;
			if (transition.get_state()[node] == past_node_value) {
				auto predictive_vector = convert_to_predictive_vector(transition.get_state();
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
