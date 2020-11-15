#ifndef OCCUPATION_TIMES_H
#define OCCUPATION_TIMES_H
#include <unordered_map>

#include "state.h"

template <class Real_t> class OccupationTimes {
private: 
	std::unordered_map<State, Real_t, StateHash> total_occupation_time;
public:
	OccupationTimes<Real_t>() {}

	void add(const State &state, Real_t time) {
		if (total_occupation_time.find(state) == total_occupation_time.end()) {
			total_occupation_time[state] = 0.0;
		}
		total_occupation_time[state] += time;
	}

	inline Real_t get_occupation_time(const State &state) const {
		return total_occupation_time.at(state);
	}
};
#endif // !OCCUPATION_TIMES_H
