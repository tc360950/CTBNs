#ifndef TRANSITION_H
#define TRANSITION_H
#include "state.h"

class Transition {
private:
	State state; 
	size_t changing_node;
	size_t new_node_state;

public:
	Transition(State state, size_t cn, size_t nv):
		state{ state },
		changing_node{ cn },
		new_node_state{ nv } {

	}
	bool operator==(const Transition &other) const
	{
		return this->state == other.state 
			&& this->changing_node == other.changing_node
			&& this->new_node_state == other.new_node_state;
	}

	const State &get_state() const {
		return this->state;
	}

	size_t get_changing_node() const {
		return this->changing_node;
	}

	size_t get_old_node_state() const {
		return state.get_node_value(changing_node);
	}

	static Transition create_from_states(const State &old_state, const State &new_state) {
		size_t changing_node = 0;
		for (size_t i = 0; i < old_state.get_size(); i++) {
			if (old_state.get_node_value(i) != new_state.get_node_value(i)) {
				changing_node = i;
				break;
			}
		}
		return Transition{ old_state, changing_node, new_state.get_node_value(changing_node) };
	}

	friend class TransitionHash;
};

class TransitionHash {
public:
	size_t operator()(const Transition& s) const
	{
		return (std::hash<std::vector<bool>>()(s.state.get_data()) + s.changing_node);
	}
};

#endif // !TRANSITION_H
