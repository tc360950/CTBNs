#ifndef STATE_H
#define STATE_H
#include <vector>

class State {
private:
	std::vector<bool> state;

public:
	State(std::vector<bool> state) : state{ state } {

	}

	State(const State &s2) { this->state = s2.state; }

	const std::vector<bool> &get_data() const {
		return state;
	}
	bool operator==(const State &other) const
	{
		return this->state == other.state;
	}

	void flip_node_value(size_t node) {
		state[node] = !state[node];
	}

	inline size_t get_size() const {
		return state.size();
	}

	inline bool get_node_value(size_t node) const {
		return state[node];
	}

	friend class StateHash;
};

class StateHash {
public:
	size_t operator()(const State& s) const
	{
		return (std::hash<std::vector<bool>>()(s.state));
	}
};

#endif // !STATE_H
