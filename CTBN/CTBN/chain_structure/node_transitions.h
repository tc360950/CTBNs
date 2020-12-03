#ifndef NODE_TRANSITIONS_H
#define NODE_TRANSITIONS_H
#include <vector>

#include "transition.h"

template<class Real_t> class NodeTransitions {
public:
	std::vector<Real_t> state_counts;
	std::vector<Real_t> time_spent_in_state;
	std::vector<Real_t> predictive_vectors;
	/**
	Pierwszy skladnik sumy, ktora oblicza gradient to suma n(c, s, s') * Z(c), w tym wektorze trzymamy te sume.
	**/
	std::vector<Real_t> sum_counts_times_predictive;
	std::vector<Real_t> predictive_vectors_times_occupation;

	size_t number_of_value_changing_jumps;
	size_t number_of_jumps_from_node_value;
	NodeTransitions<Real_t>() {
	}

	void add(const std::vector<Real_t> &predictive, const Real_t time, const Real_t count) {
		state_counts.push_back(count);
		time_spent_in_state.push_back(time);
		for (size_t i = 0; i < predictive.size(); i++) {
			predictive_vectors.push_back(predictive[i]);
		}
	}

	void end_add(const size_t parameters_size) {
		sum_counts_times_predictive.resize(parameters_size);
		predictive_vectors_times_occupation.resize(predictive_vectors.size());
		for (auto &el : sum_counts_times_predictive) {
			el = 0.0;
		}
		size_t start = 0;
		for (size_t i = 0; i < time_spent_in_state.size(); i++) {
			for (size_t j = 0; j < parameters_size; j++) {
				predictive_vectors_times_occupation[start + j] = predictive_vectors[start + j] * time_spent_in_state[i];
				sum_counts_times_predictive[j] -= state_counts[i] * predictive_vectors[start + j];
			}
			start += parameters_size;
		}
	}

	void gather_jump_information(const size_t node, const size_t node_value, const std::vector<std::pair<State, Real_t>> &skeleton) {
		for (size_t i = 0; i < skeleton.size() - 1; i++) {
			if (skeleton[i].first.get_node_value(node) == node_value) {
				number_of_jumps_from_node_value++;
				if (skeleton[i].first.get_node_value(node) != skeleton[i + 1].first.get_node_value(node)) {
					number_of_value_changing_jumps++;
				}
			}
		}
	}
};

#endif // !NODE_TRANSITIONS_H