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

	std::string vector_to_string(const std::vector<Real_t> &vec) {
		std::string vector;
		for (auto v : vec) {
			vector.append(std::to_string(v));
			vector.append(" ");
		}
		return vector;
	}
public:

	void random_test(long model_seed) {
		logTest<LikelihoodTest>("Starting likelihood test");
		EmptyModel<Real_t> model{ 20, model_seed };
		auto model_data = model.sample_chain_and_skeleton(T_MAX);
		LikelihoodCalculator<Real_t> calculator{ model_data.first.transition_repository, T_MAX };
		LikelihoodTester<Real_t> tester{ NUMBER_OF_NODES, T_MAX };
		auto beta = sample_beta();
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

		for (size_t i = 0; i < 2 * NUMBER_OF_NODES; i++) {
			auto node = i / 2;
			auto value = i % 2;
			std::vector<Real_t> likelihood_gradient_1;
			likelihood_gradient_1.resize(NUMBER_OF_NODES);
			calculator.calculate_likelihood_gradient(beta, node, value, likelihood_gradient_1);
			auto likelihood_gradient_2 = tester.calculate_gradient(beta, model_data.second, node, value);
			bool ok = true;
			for (size_t c = 0; c < likelihood_gradient_2.size(); c++) {
				if (std::abs(likelihood_gradient_1[c] - likelihood_gradient_2[c]) >= TOLERANCE) {
					logTest<LikelihoodTest>("Likelihood comparison fail with: ", vector_to_string(likelihood_gradient_1), " ", vector_to_string(likelihood_gradient_2));
					ok = false;
				}
			}
			if (ok) {
				logTest<LikelihoodTest>("Likelihood comparison ok with: ", vector_to_string(likelihood_gradient_1), " ||| ", vector_to_string(likelihood_gradient_2));
			}
		}
	}
};
#endif // !LIKELIHOOD_TEST_H
