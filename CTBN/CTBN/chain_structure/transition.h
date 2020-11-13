#ifndef TRANSITION_H
#define TRANSITION_H
#include "state.h"
//TODO getters, new State formation, hash function
class Transition {
private:
	State state; 
	size_t changing_node;
	bool new_node_state;

public:
	bool operator==(const Transition &other) const
	{
		return this->state == other->state 
			&& this->changing_node == other->changing_node
			&& this->new_node_state == other.new_node_state;
	}

	const State &get_state() const {
		return this->state;
	}
};

#endif // !TRANSITION_H
