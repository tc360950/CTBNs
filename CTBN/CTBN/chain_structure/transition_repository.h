#ifndef TRANSITION_REPOSITORY_H
#define TRANSITION_REPOSITORY_H
#include <vector>

#include "node_transitions.h"
#include "occupation_times.h"

template <class Real_t> class TransitionRepository {
private:
	std::vector<NodeTransitions<Real_t>> node_transitions;
	const size_t number_of_nodes;
	const size_t parameters_size;

public:
	TransitionRepository<Real_t>(std::vector<NodeTransitions<Real_t>> nt, size_t nn, size_t ps):
		node_transitions{nt},
		number_of_nodes{nn},
		parameters_size{ ps }{}

	const NodeTransitions<Real_t> &fetch_node_transitions(const size_t node, const size_t past_node_value) const {
		return node_transitions[2*node + past_node_value];
	}
	const size_t get_parameters_size() const {
		return parameters_size;
	}
};

#endif // !TRANSITION_REPOSITORY_H
