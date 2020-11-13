#ifndef NODE_TRANSITIONS_H
#define NODE_TRANSITIONS_H
#include <unordered_map>

#include "transition.h"

template<class Real_t> class NodeTransitions {
private:
	std::unordered_map<Transition, Real_t> transition_counts;
	size_t node;

};

#endif // !NODE_TRANSITIONS_H
