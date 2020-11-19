#ifndef LIKELIHOOD_TEST_H
#define LIKELIHOOD_TEST_H
#include <random>

#include "../../models/empty_model.h"
#include "likelihood_tester.h"
#include "../likelihood_calculator.h"
#include "../../utils/logger.h"

template <class Real_t> class LikelihoodTest {
private:
	const size_t NUMBER_OF_NODES = 20;
	const long SEED = 213141;
	const Real_t T_MAX = 50.0;
	const Real_t TOLERANCE = 0.0001;
	std::mt19937 generator{SEED};
	std::vector<Real_t> sample_beta() {
		std::vector<Real_t> result;
		std::uniform_real_distribution<Real_t> distribution(0.0, 1.0);
		for (size_t i = 0; i < NUMBER_OF_NODES; i++) {
			result.push_back(distribution(generator));
		}
		return result;
	}
public:

	void random_test(long model_seed) {
		logTest<LikelihoodTest>("Starting likelihood test");
		EmptyModel<Real_t> model{ 20, model_seed };
		auto model_data = model.sample_chain_and_skeleton(T_MAX);
		LikelihoodCalculator<Real_t> calculator{ model_data.first.transition_repository, T_MAX };
		LikelihoodTester<Real_t> tester{ NUMBER_OF_NODES, T_MAX };
		auto beta = sample_beta();
		logTest<LikelihoodTest>("DUPA");
		for (size_t i = 0; i < 2 * NUMBER_OF_NODES; i++) {
			auto node = i / 2;
			auto value = i % 2;
			auto likelihood_1 = calculator.calculate_likelihood(beta, node, value);
			auto likelihood_2 = tester.calculate_likelihood(beta, model_data.second, node, value);
			if (std::abs(likelihood_1 - likelihood_2) >= TOLERANCE) {
				logTest<LikelihoodTest>("Likelihood comparison fail with: ", likelihood_1, " ", likelihood_2);
			}
			else {
				logTest<LikelihoodTest>("Likelihood comparison ok with: ", likelihood_1, " ", likelihood_2);
			}
		}
	}
};
#endif // !LIKELIHOOD_TEST_H
