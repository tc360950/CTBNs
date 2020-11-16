#ifndef NODE_TRANSITIONS_H
#define NODE_TRANSITIONS_H
#include <vector>

#include "transition.h"

template<class Real_t> class NodeTransitions {
public:
	std::vector<Real_t> state_counts;
	std::vector<Real_t> time_spent_in_state;
	std::vector<Real_t> predictive_vectors;

	NodeTransitions<Real_t>() {
	}

	void add(const std::vector<Real_t> &predictive, const Real_t time, const Real_t count) {
		state_counts.push_back(count);
		time_spent_in_state.push_back(time);
		for (auto el : predictive) {
			predictive_vectors.push_back(el);
		}
	}

};

#endif // !NODE_TRANSITIONS_H
