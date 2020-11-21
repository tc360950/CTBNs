#ifndef EMPTY_MODEL_H
#define EMPTY_MODEL_H
#include <random>
#include <utility>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../chain_structure/transition_repository.h"
#include "model_data.h"
#include "../utils/logger.h"

template <class Real_t> class EmptyModel {
public:
	std::mt19937 generator;
	std::vector<bool> preferences;

	State sample_state(const State &current_state, const std::vector<Real_t> &intensities) {
		State new_state{ current_state };
		std::discrete_distribution<size_t> distribution(intensities.begin(), intensities.end());
		auto node = distribution(generator);
		new_state.flip_node_value(node);
		return new_state;
	}

	std::vector<std::pair<State, Real_t>> simulate(Real_t t_max, const State &starting_state) {
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
		skeleton.pop_back();
		if (TEST) {
			for (size_t i = 0; i < skeleton.size() - 1; i++) {
				if (skeleton[i].second >= skeleton[i + 1].second) {
					logTest<EmptyModel>("EmptyModel: Jump times are not increasing!");
				}
				size_t changing_count = 0;
				for (size_t j = 0; j < skeleton[i].first.get_size(); j++) {
					if (skeleton[i].first.get_node_value(j) != skeleton[i + 1].first.get_node_value(j)) {
						changing_count++;
					}
				}
				if (changing_count != 1) {
					logTest<EmptyModel>("More or less than one node has been changed in a jump!");
				}
			}
			logTest<EmptyModel>("Skeleton is ok!");
		}
		return skeleton;
	}

	Real_t get_intensity(size_t node, bool node_value, const State &state) const {
			return 5.0;
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

	std::vector<Real_t> state_to_predictive_vector(const State &state, const size_t node) const {
		std::vector<Real_t> result;
		result.push_back(1.0);
		for (size_t i = 0; i < state.get_size(); i++) {
			if (i != node) {
				result.push_back(state.get_data()[i]);
			}
		}
		return result;
	}

	Real_t extract_transition_count(const State &state, std::unordered_map<State, std::unordered_set<Transition, TransitionHash>, StateHash> &state_to_transition,
		std::unordered_map<Transition, Real_t, TransitionHash> &transition_to_count, const size_t node, const size_t node_past_value) const {
		for (auto &transition : state_to_transition[state]) {
			if (transition.get_changing_node() == node && transition.get_old_node_state() == node_past_value) {
				return transition_to_count[transition];
			}
		}
		return 0.0;
	}

	TransitionRepository<Real_t> convert_skeleton_to_transition_repository(const std::vector<std::pair<State, Real_t>> &skeleton, const Real_t t_max) const {
		OccupationTimes<Real_t> occupation_times = extract_occupation_times(skeleton, t_max);
		std::vector<NodeTransitions<Real_t>> node_transitions;
		std::unordered_map<Transition, Real_t, TransitionHash> transition_to_count;
		std::unordered_map<State, std::unordered_set<Transition, TransitionHash>, StateHash> state_to_transition;
		auto all_states = occupation_times.get_states();

		//PREPROCESSING
		for (size_t i = 0; i < 2 * preferences.size(); i++) {
			node_transitions.push_back(NodeTransitions<Real_t>());
		}
		for (auto state : all_states) {
			state_to_transition[state] = std::unordered_set<Transition, TransitionHash>();
		}

		for (size_t i = 1; i < skeleton.size(); i++) {
			auto transition = Transition::create_from_states(skeleton[i - 1].first, skeleton[i].first);
			state_to_transition[skeleton[i - 1].first].insert(transition);
			if (transition_to_count.find(transition) == transition_to_count.end()) {
				transition_to_count[transition] = 0.0;
			}
			transition_to_count[transition]++;
		}
		//END PREPROCESSING
		if (TEST) {
			if (all_states.size() > skeleton.size()) {
				logTest<EmptyModel>("There are more states in <all_states> than in the skeleton!");
			}
			std::unordered_set<State, StateHash> all_states_set(all_states.begin(), all_states.end());
			if (all_states.size() != all_states_set.size()) {
				logTest<EmptyModel>("There are duplicate states in <all_states>!");
			}
			Real_t total_transition_count = 0.0;
			for (auto &tran : transition_to_count) {
				total_transition_count += tran.second;
			}
			if (std::abs(total_transition_count - skeleton.size() + 1.0) > 0.5) {
				logTest<EmptyModel>("There are more transition counts than states!");
			}
			size_t total_transitions = 0;
			size_t zero_transition_per_state = 0;
			for (auto &el : state_to_transition) {
				total_transitions += el.second.size();
				if (el.second.size() == 0) {
					zero_transition_per_state++;
				}
			}
			if (total_transitions != transition_to_count.size()) {
				logTest<EmptyModel>("Transition counts do not match!");
			}
			if (zero_transition_per_state > 1) {
				logTest<EmptyModel>("There are more than one state with zero transitions!");
			}
		}

		for (size_t i = 0; i < 2 * preferences.size(); i++) {
			const size_t node = i / 2;
			const size_t past_node_value = i % 2;
			for (auto & state : all_states) {
				if (state.get_node_value(node) == past_node_value) {
					auto time = occupation_times.get_occupation_time(state);
					auto predictive_vector = state_to_predictive_vector(state, node);
					auto transition_count = extract_transition_count(state, state_to_transition, transition_to_count, node, past_node_value);
					node_transitions[i].add(predictive_vector, time, transition_count);
				}
			}
		}

		for (auto &nt : node_transitions) {
			nt.end_add(preferences.size());
		}
		return TransitionRepository<Real_t>{node_transitions, preferences.size(), preferences.size()};
	}

	std::vector<std::vector<bool>> generate_dependence_structure() const {
		std::vector<std::vector<bool>> dependence{ preferences.size() };
		for (auto &vec : dependence) {
			for (size_t i = 0; i < preferences.size(); i++) {
				vec.push_back(false);
			}
		}
		return dependence;
	}

	State simulate_starting_state() {
		std::vector<bool> state;
		for (size_t i = 0; i < preferences.size(); i++) {
			state.push_back(random_bit());
		}
		return State(state);
	}

	bool random_bit() {
		std::uniform_int_distribution<int> distrib(0, 1);
		auto bit = distrib(generator);
		return  bit == 1 ? true : false;
	}

	//Dor test purposes only!
	std::pair<ModelData<Real_t>, std::vector<std::pair<State, Real_t>>> sample_chain_and_skeleton(Real_t t_max) {
		auto starting_state = simulate_starting_state();
		auto skeleton = simulate(t_max, starting_state);
		log("Simulated skeleton for list model of size ", skeleton.size());
		auto transitions = convert_skeleton_to_transition_repository(skeleton, t_max);
		log("Converted skeleton to transition repository");
		auto dependence_structure = generate_dependence_structure();
		return std::make_pair(ModelData<Real_t>(transitions, skeleton.size(), dependence_structure), skeleton);
	}

public:
	EmptyModel<Real_t>(size_t number_of_nodes, long seed) :
		generator{ seed } {
		preferences.resize(number_of_nodes);
		for (size_t i = 0; i < number_of_nodes; i++) {
			preferences[i] = random_bit();
		}
	}

	ModelData<Real_t> sample_chain(Real_t t_max) {
		auto starting_state = simulate_starting_state();
		auto skeleton = simulate(t_max, starting_state);
		log("Simulated skeleton for list model of size ", skeleton.size());
		auto transitions = convert_skeleton_to_transition_repository(skeleton, t_max);
		log("Converted skeleton to transition repository");
		auto dependence_structure = generate_dependence_structure();
		return ModelData<Real_t>(transitions, skeleton.size(), dependence_structure);
	}
};

#endif // !EMPTY_MODEL_H
