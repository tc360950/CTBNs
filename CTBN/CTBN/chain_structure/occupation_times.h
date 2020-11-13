#ifndef OCCUPATION_TIMES_H
#define OCCUPATION_TIMES_H
#include <unordered_map>

#include "state.h"

template <class Real_t> class OccupationTimes {
private: 
	std::unordered_map<State, Real_t> total_occupation_time;
public:
	Real_t get_occupation_time(const State &state) const {
		return total_occupation_time[state];
	}
};
#endif // !OCCUPATION_TIMES_H
