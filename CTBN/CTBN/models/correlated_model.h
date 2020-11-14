#ifndef CORRELATED_MODEL_H
#define CORRELATED_MODEL_H
#include <random>
#include <pair>

#include "../chain_structure/transition_repository.h"

template <class Real_t> class CorrelatedModel {
private:
	std::mt19937 generator;
	//tylko pierwsze 5 nodeow
	std::vector<std::pair<size_t, size_t>> parents;
	std::vector<bool> preferences;

	State sample_state(const State &current_state, const std::vector<Real_t> &intensities) const {
		State new_state{ current_state };
		std::discrete_distribution<Real_t> distribution(intensities);
		auto node = distribution(generator);
		new_state.flip_node_value(node);
		return new_state;
	}

	std::vector<std::pair<State, Real_t>> simulate(Real_t t_max, const State &starting_state) const {
		//stan + czas skoku do stanu 
		std::vector<std::pair<State, Real_t>> skeleton;
		skeleton.push_back(std::make_pair(starting_state, 0.0));
		while (skeleton.back().second <= t_max) {
			auto current_state = skeleton.back().first;
			auto intensities = create_intensity_weights(current_state);
			Real_t state_intensity = std::accumulate(intensities.begin(), intensities.end(), 0.0);
			std::exponential_distribution<Real_t> exponential(state_intensity);
			auto time = exponential(generator) + skeleton.back().second;
			auto new_state = sample_state(current_state, intensities);
			skeleton.push_back(std::make_pair(new_state, time));
		}
		skeleton.erase(skeleton.size());
		return skeleton;
	}

	Real_t get_intensity(size_t node, bool node_value, const State &state) const {
		if (node > 4) {
			return 5.0;
		}
		else {
			auto parent_product = state.get_node_value(parents[node].first) && state.get_node_value(parents[node].second);
			if (state.get_node_value(node) == preferences[node]) {
				return parent_product ? 9.0 : 1.0;
			}
			else {
				return parent_product ? 1.0 : 0.0;
			}
		}
	}

	std::vector<Real_t> create_intensity_weights(const State &state) const {
		std::vector<Real_t> result(state.get_size());
		for (size_t i = 0; i < state.get_size(); i++) {
			result[i] = get_intensity(i, state.get_node_value(i), state);
		}
		return result;
	}
	OccupationTimes<Real_t> extract_occupation_times(const std::vector<std::pair<State, Real_t>> &skeleton, const Real_t t_max) const {
		OccupationTimes<Real_t> occupation_times;
		for (size_t i = 0; i < skeleton.size() - 1; i++) {
			occupation_times.add(skeleton[i].first, skeleton[i + 1].second - skeleton[i].second);
		}
		occupation_times.add(skeleton.back().first, t_max - skeleton.back().second);
		return occupation_times;
	}

	TransitionRepository<Real_t> convert_skeleton_to_transition_repository(const std::vector<std::pair<State, Real_t>> &skeleton, const Real_t t_max) const {
		OccupationTimes<Real_t> occupation_times extract_occupation_times(skeleton, t_max);
		std::vector<NodeTransitions<Real_t>> node_transitions;
		for (size_t i = 0; i < preferences.size(); i++) {
			node_transitions.push_back(NodeTransitions<Real_t>(i));
		}
		for (size_t i = 1; i < skeleton.size(); i++) {
			auto transition = Transition::create_from_states(skeleton[i - 1].first, skeleton[i].first);
			node_transitions[transition.get_changing_node()].add(transition);
		}
		return TransitionRepository<Real_t>(occupation_times, node_transitions);
	}

public:
	CorrelatedModel<Real_t>(size_t number_of_nodes, long seed) :
		generator{ seed },
		preferences{ number_of_nodes } {
		std::uniform_int_distribution<int> distrib(0, 1);
		for (size_t i = 0; i < number_of_nodes; i++) {
			auto bit = distrib(generator);
			preferences[i] = bit == 1 ? true : false;
		}
		std::uniform_int_distribution<size_t> distribution(0, 4);
		for (size_t i = 0; i < 5; i++) {
			auto parent_1 = distribution(generator);
			auto parent_2 = distribution(generator);
			while (parent_1 == parent_2 || parent_1 == i || parent_2 == i) {
				auto parent_1 = distribution(generator);
				auto parent_2 = distribution(generator);
			}
			parents.push_back(std::make_pair(parent_1, parent_2));
		}
	}

	//number of jumps is returned 
	std::pair<TransitionRepository<Real_t>, Real_t> sample_chain(Real_t t_max, const State &starting_state) const {
		auto skeleton = simulate(t_max, starting_state);
		return convert_skeleton_to_transition_repository(skeleton, t_max);
	}
};
#endif // !CORRELATED_MODEL_H
