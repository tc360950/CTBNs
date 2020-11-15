#ifndef MODEL_DATA_H
#define MODEL_DATA_H

#include "../chain_structure/transition_repository.h"

template <class Real_t> class ModelData {
public:
	const TransitionRepository<Real_t> transition_repository;
	const Real_t number_of_jumps;
	const std::vector<std::vector<bool>> dependence_structure;

	ModelData<Real_t>(TransitionRepository<Real_t> tr, Real_t nj, std::vector<std::vector<bool>> ds):
		transition_repository {tr},
		number_of_jumps{ nj },
		dependence_structure{ ds } {

	}
};

#endif // !MODEL_DATA_H
