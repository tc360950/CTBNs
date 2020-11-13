#ifndef STATE_H
#define STATE_H
#include <vector>

class State {
private:
	std::vector<bool> state;

public:
	bool operator==(const State &other) const
	{
		return this->state == other->state;
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
		return (std::hash<std::vector<bool>>()(s->state));
	}
};

#endif // !STATE_H
