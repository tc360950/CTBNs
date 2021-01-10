#include <random>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <ctime>

#include "CTBN.h"
#include "bob_dylan.h"
#include "models/list_model.h"
#include "models/correlated_model.h"
#include "summary_statistics/statistics_factory.h"
#include "models/empty_model.h"
#include "models/correlated_model_no_interactions.h"
#include "likelihood_calculator/tests/likelihood_test.h"
#include "solvers/tests/admm_solver_test.h"
#include "utils/logger.h"
#include "utils/random_tree_generator.h"
#include "models/det_binary_tree.h"
#include "utils/parameters.h"

template <class Model> void simulate_memory(double t_max, size_t no_of_nodes, const std::string NAME) {
	const size_t TRIES = 100;
	std::vector<Statistics<double>> results;
	results.resize(TRIES);
	std::srand(std::time(nullptr));

	for (size_t i = 0; i < TRIES; i++) {
		auto seed = std::rand();
		BobDylan<double, Model> bob;
		StatisticsFactory<double, Model> factory;
		auto result = bob.simulate(no_of_nodes, seed, t_max);
		auto stats = factory.convert(no_of_nodes, t_max, result.model_data, result);
		results[i] = stats;
	}
	for (size_t i = 1; i < TRIES; i++) {
		results[0].add(results[i]);
	}
	logTest<Model>(NAME, " Results: \n", "FDR: ", results[0].FDR / (TRIES), "\nMD: ", results[0].MD/ ( TRIES), "\nPOWER: ", results[0].power / (TRIES));
}



void test(long seed) {
	//Works only for MEMORY_HUNGRY 
	{
		ADMMSolverTest<double, ListModel<double>> admm_test;
		LikelihoodTest<double, ListModel<double>> tester;
		for (size_t i = 0; i < 10; i++) {
			admm_test.test(seed, 0.1);
			tester.random_test(seed);
			seed++;
		}
	}
	{
		ADMMSolverTest<double, CorrelatedModelNoInteractions<double>> admm_test;
		LikelihoodTest<double, CorrelatedModelNoInteractions<double>> tester;
		for (size_t i = 0; i < 10; i++) {
			admm_test.test(seed, 0.1);
			tester.random_test(seed);
			seed++;
		}
	}
}

int main(int argc, char **argv)
{	
	if (argc != 2) {
		return -1;
	}
	long seed = std::stol(argv[1]);

		simulate_memory<ListModel<double>>(10, 20, "10_20");
		simulate_memory<ListModel<double>>(50, 20, "50_20");
		simulate_memory<CorrelatedModelNoInteractions<double>>(10, 20, "10_20");
		simulate_memory<CorrelatedModelNoInteractions<double>>(50, 20, "50_20");
		simulate_memory<CorrelatedModel<double>>(10, 20, "10_20");
		simulate_memory<CorrelatedModel<double>>(50, 20, "50_20");
		simulate_memory<BinaryTree<double>>(10, 20, "10_20");
		simulate_memory<BinaryTree<double>>(50, 20, "50_20");

		simulate_memory<ListModel<double>>(10, 50, "10_50");
		simulate_memory<ListModel<double>>(50, 50, "50_50");
		simulate_memory<CorrelatedModelNoInteractions<double>>(10, 50, "10_50");
		simulate_memory<CorrelatedModelNoInteractions<double>>(50, 50, "50_50");
		simulate_memory<CorrelatedModel<double>>(10, 50, "10_50");
		simulate_memory<CorrelatedModel<double>>(50, 50, "50_50");
		simulate_memory<BinaryTree<double>>(10, 50, "10_50");
		simulate_memory<BinaryTree<double>>(50, 50, "50_50");

	return 0;
}
