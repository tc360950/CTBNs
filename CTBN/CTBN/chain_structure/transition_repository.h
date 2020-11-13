#ifndef TRANSITION_REPOSITORY_H
#define TRANSITION_REPOSITORY_H
#include <vector>

#include "node_transitions.h"
#include "occupation_times.h"

template <class Real_t> class TransitionRepository {
private:
	//pod indeksem i tranzycje i-tego wierzchola
	std::vector<NodeTransitions<Real_t>> node_transitions;
	OccupationTimes<Real_t> occupation_times;

public:
	const NodeTransitions<Real_t> &fetch_node_transitions(size_t node) const {
		return node_transitions[node];
	}

	Real_t get_occupation_time(const State &state) const {
		return occupation_times.get_occupation_time(state);
	}
};

#endif // !TRANSITION_REPOSITORY_H
