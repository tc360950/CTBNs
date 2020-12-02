// CTBN.cpp: definiuje punkt wejścia dla aplikacji.
//
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
#include "models/binary_tree_model.h"

template <class Model> void simulate(double t_max, size_t no_of_nodes) {
	const size_t TRIES = 1;
	std::vector<Statistics<double>> results;
	results.resize(TRIES);
	std::vector<std::thread> threads;
	std::srand(std::time(nullptr));

	for (size_t i = 0; i < TRIES; i++) {
		auto seed = std::rand();
		threads.emplace_back([t_max, no_of_nodes, seed, i, &results ] {
			BobDylan<double, Model> bob;
			StatisticsFactory<double, Model> factory;
			auto result = bob.simulate(no_of_nodes, seed, t_max);
			auto stats = factory.convert(no_of_nodes, t_max, result.model_data, result);
			results[i] = stats;
		});
	}
	for (auto &th : threads)
	{
		th.join();
	}
	for (size_t i = 1; i < TRIES; i++) {
		results[0].add(results[i]);
	}
	logTest<Model>("Results: \n", "FDR: ", results[0].FDR / TRIES, "\nMD: ", results[0].MD/ TRIES, "\nPOWER: ", results[0].power / TRIES);
}

int main(int argc, char **argv)
{	
	long seed = 12312;//std::stol(argv[1]);
	RandomTree<double> randomTree(12312);
	auto tree = randomTree.generate_random_binary_tree(5);
	simulate<TreeModelNoInteractions<double>>(10, 20);
	simulate<TreeModelNoInteractions<double>>(50, 20);
	/*
	simulate<ListModel<double>>(10, 20);
	simulate<ListModel<double>>(50, 20);
	simulate<ListModel<double>>(50, 50);
	simulate<ListModel<double>>(10, 50);

	simulate<CorrelatedModelNoInteractions<double>>(10, 20);
	simulate<CorrelatedModelNoInteractions<double>>(50, 20);
	*/
	
	
	/*ADMMSolverTest<double> admm_test;
	//admm_test.test(1243453, 0.1);

	LikelihoodTest<double> tester;
	//tester.random_test(12188845);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	BobDylan<double, CorrelatedModelNoInteractions<double>> bob;
	auto result = bob.simulate(20, seed, 50);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[µs]" << std::endl;
	StatisticsFactory<double,CorrelatedModelNoInteractions<double>> factory;
	auto stats = factory.convert(20, 50, result.model_data, result);
    std::cout << stats.power << "\n";
    std::cout << stats.FDR << "\n";
    std::cout << stats.MD << "\n";
	while (true) {}*/
	return 0;
}
