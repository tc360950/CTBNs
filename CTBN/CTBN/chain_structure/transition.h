#ifndef TRANSITION_H
#define TRANSITION_H
#include "state.h"
//TODO getters, new State formation, hash function, == operator
class Transition {
private:
	State state; 
	size_t changing_node;
	bool new_node_state;
};

#endif // !TRANSITION_H
