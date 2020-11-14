#ifndef NODE_TRANSITIONS_H
#define NODE_TRANSITIONS_H
#include <unordered_map>

#include "transition.h"

template<class Real_t> class NodeTransitions {
private:
	std::unordered_map<Transition, Real_t, TransitionHash> transition_counts;
	size_t node;

public:
	NodeTransitions<Real_t>(size_t node) : node{ node } {
	}

	void add(const Transition &transition) {
		if (transition_counts.find(transition) == transition.end()) {
			transition_counts[transition] = 0.0;
		}
		transition_counts[transition] += 1.0;
	}

	const std::unordered_map<Transition, Real_t> &get_transition_counts() const {
		return this->transition_counts;
	}

};

#endif // !NODE_TRANSITIONS_H
