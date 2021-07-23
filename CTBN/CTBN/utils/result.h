#ifndef SIMULATION_RESULT_H
#define SIMULATION_RESULT_H
#include <vector>

#include "../models/model_data.h"

template <class Real_t> class Result {
public:
	std::vector<Real_t> beta;
	size_t node;
	size_t node_past_value;
	Result<Real_t>(std::vector<Real_t> b, size_t n, size_t nv) :
		beta{ b },
		node{ n },
		node_past_value{ nv } {

	}
	Result<Real_t>() {}
};

template <class Real_t> struct SimulationResult {
	std::vector<Result<Real_t>> inference_results;
	ModelData<Real_t> model_data;


};
#endif // !SIMULATION_RESULT_H
