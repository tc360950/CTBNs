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
#include "models/det_binary_tree.h"
#include "utils/parameters.h"

template <class Model> void simulate(double t_max, size_t no_of_nodes, const std::string NAME) {
	const size_t TRIES = 50;
	std::vector<Statistics<double>> results;
	results.resize(2*TRIES);
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
	threads.clear();
	for (size_t i = 0; i < TRIES; i++) {
		auto seed = std::rand();
		threads.emplace_back([t_max, no_of_nodes, seed, i, &results, TRIES] {
			BobDylan<double, Model> bob;
			StatisticsFactory<double, Model> factory;
			auto result = bob.simulate(no_of_nodes, seed, t_max);
			auto stats = factory.convert(no_of_nodes, t_max, result.model_data, result);
			results[i + TRIES] = stats;
		});
	}
	for (auto &th : threads)
	{
		th.join();
	}


	for (size_t i = 1; i < 2 * TRIES; i++) {
		results[0].add(results[i]);
	}
	logTest<Model>(NAME, " Results: \n", "FDR: ", results[0].FDR / (2* TRIES), "\nMD: ", results[0].MD/ (2* TRIES), "\nPOWER: ", results[0].power / (2* TRIES));
}


template <class Model> void simulate_memory(double t_max, size_t no_of_nodes, const std::string NAME) {
	const size_t TRIES = 25;
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
	logTest<Model>(NAME, " Results: \n", "FDR: ", results[0].FDR / (TRIES), "\nMD: ", results[0].MD/ ( TRIES), "\nPOWER: ", results[0].power / ( TRIES));
}



void test(long seed) {
	{ADMMSolverTest<double, ListModel<double>> admm_test;

	LikelihoodTest<double, ListModel<double>> tester;
	for (size_t i = 0; i < 10; i++) {
		admm_test.test(seed, 0.1);
		tester.random_test(seed);
		seed++;
	}}
	{ADMMSolverTest<double, CorrelatedModelNoInteractions<double>> admm_test;

	LikelihoodTest<double, CorrelatedModelNoInteractions<double>> tester;
	for (size_t i = 0; i < 10; i++) {
		admm_test.test(seed, 0.1);
		tester.random_test(seed);
		seed++;
	}}
}

int main(int argc, char **argv)
{	
	if (argc != 3) {
		return -1;
	}
	long seed = std::stol(argv[1]);
	N_DEFINITION = std::stoi(argv[2]);
	//long seed = 21435423;
	//N_DEFINITION = 0;
	//test(seed);
	//test(seed);


	{N_DEFINITION = 0;
	//	simulate<ListModel<double>>(10, 20, "10_20");
		//simulate<ListModel<double>>(50, 20, "50_20");
		//simulate<CorrelatedModelNoInteractions<double>>(10, 20, "10_20");
		//simulate<CorrelatedModelNoInteractions<double>>(50, 20, "50_20");
		//simulate<CorrelatedModel<double>>(10, 20, "10_20");
		//simulate<CorrelatedModel<double>>(50, 20, "50_20");
		//simulate<BinaryTree<double>>(10, 20, "10_20");
		//simulate<BinaryTree<double>>(50, 20, "50_20");

		//simulate<ListModel<double>>(10, 50, "10_50");
		//simulate<ListModel<double>>(50, 50, "50_50");
		//simulate<CorrelatedModelNoInteractions<double>>(10, 50, "10_50");
		//simulate<CorrelatedModelNoInteractions<double>>(50, 50, "50_50");
		simulate_memory<CorrelatedModel<double>>(10, 50, "10_50");
		simulate_memory<CorrelatedModel<double>>(50, 50, "50_50");
		//simulate<BinaryTree<double>>(10, 50, "10_50");
		//simulate<BinaryTree<double>>(50, 50, "50_50");
return 0;
	}
	std::cout << "STOP\n\n\n\n";
	{N_DEFINITION = 3;
	simulate<ListModel<double>>(10, 20, "10_20");
	simulate<ListModel<double>>(50, 20, "50_20");
	simulate<CorrelatedModelNoInteractions<double>>(10, 20, "10_20");
	simulate<CorrelatedModelNoInteractions<double>>(50, 20, "50_20");
	simulate<CorrelatedModel<double>>(10, 20, "10_20");
	simulate<CorrelatedModel<double>>(50, 20, "50_20");
	simulate<BinaryTree<double>>(10, 20, "10_20");
	simulate<BinaryTree<double>>(50, 20, "50_20");

	simulate<ListModel<double>>(10, 50, "10_50");
	simulate<ListModel<double>>(50, 50, "50_50");
	simulate<CorrelatedModelNoInteractions<double>>(10, 50, "10_50");
	simulate<CorrelatedModelNoInteractions<double>>(50, 50, "50_50");
	simulate<CorrelatedModel<double>>(10, 50, "10_50");
	simulate<CorrelatedModel<double>>(50, 50, "50_50");
	simulate<BinaryTree<double>>(10, 50, "10_50");
	simulate<BinaryTree<double>>(50, 50, "50_50");
	}
	std::cout << "STOP\n\n\n\n";

/*	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
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
