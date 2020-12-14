#ifndef LIST_MODEL_H
#define LIST_MODEL_H
#include <random>
#include <utility>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../chain_structure/transition_repository.h"
#include "model_data.h"
#include "../utils/logger.h"

template <class Real_t> class ListModel {
private:
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
		return skeleton;
	}

	Real_t get_intensity(size_t node, bool node_value, const State &state) const {
		if (node == 0) {
			return 5.0;
		}
		else if (state.get_node_value(node - 1) == preferences[node]) {
			return node_value ? 9.0 : 1.0;
		}
		else {
			return node_value ? 1.0 : 9.0;
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

	//TODO mozna zoptymalizowac przy pomocy move semantics, jezeli bedzie za wole
	TransitionRepository<Real_t> convert_skeleton_to_transition_repository(const std::vector<std::pair<State, Real_t>> &skeleton, const Real_t t_max) const {
		OccupationTimes<Real_t> occupation_times = extract_occupation_times(skeleton, t_max);
		std::vector<NodeTransitions<Real_t>> node_transitions;
		std::unordered_map<Transition, Real_t, TransitionHash> transition_to_count;
		std::unordered_map<State, std::unordered_set<Transition, TransitionHash>, StateHash> state_to_transition;

		for (size_t i = 0; i < 2 * preferences.size(); i++) {
			node_transitions.push_back(NodeTransitions<Real_t>());
		}
		auto all_states = occupation_times.get_states();
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

		for (size_t i = 0; i < 2 * preferences.size(); i++) {
			const size_t node = i / 2;
			const size_t past_node_value = i % 2;
			for (auto & state : all_states) {
				if (state.get_node_value(node) == past_node_value) {
					auto time = occupation_times.get_occupation_time(state);
					auto predictive_vector = state_to_predictive_vector(state, node);
					auto transition_count = extract_transition_count(state, state_to_transition, transition_to_count, node, past_node_value);
					node_transitions[2 * node + past_node_value].add(predictive_vector, time, transition_count);
				}
			}
		}

		for (auto &nt : node_transitions) {
			nt.end_add(preferences.size());
		}
		for (size_t i = 0; i < 2 * preferences.size(); i++) {
			node_transitions[i].gather_jump_information(i / 2, i % 2, skeleton);
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
		for (size_t i = 0; i < preferences.size() - 1; i++) {
			dependence[i][i + 1] = true;
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
public:
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
	ListModel<Real_t>(size_t number_of_nodes, long seed): 
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

#endif // !LIST_MODEL_H
