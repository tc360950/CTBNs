// CTBN.cpp: definiuje punkt wejścia dla aplikacji.
//
#include <random>
#include <chrono>


#include "CTBN.h"
#include "bob_dylan.h"
#include "models/list_model.h"
#include "models/correlated_model.h"
#include "summary_statistics/statistics_factory.h"
#include "models/empty_model.h"
#include "models/correlated_model_no_interactions.h"
#include "likelihood_calculator/tests/likelihood_test.h"
#include "solvers/tests/admm_solver_test.h"
int main(int argc, char **argv)
{	
    long seed = std::stol(argv[1]);
	ADMMSolverTest<double> admm_test;
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
	while (true) {}
	return 0;
}
