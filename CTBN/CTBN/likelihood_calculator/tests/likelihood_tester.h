#ifndef LIKELIHOOD_TESTER_H
#define LIKELIHOOD_TESTER_H
#include <utility>
#include <vector>

#include "../../chain_structure/state.h"
//Only for empty model
template <class Real_t> class LikelihoodTester {

private:
	const size_t number_of_nodes;
	const Real_t t_max;

	Real_t get_intensity(const std::vector<Real_t> &beta, const State &state, const size_t node) {
		Real_t result = beta[0];
		for (size_t i = 0; i < number_of_nodes; i++) {
			if (i != node) {
				if (state.get_node_value(i) == 1) {
					const size_t index = i < node ? i : i - 1;
					result += beta[index + 1];
				}
			}
		}
		return result;
	}

	std::vector<Real_t> get_predictive_vector(const State &state, const size_t node) {
		std::vector<Real_t> result;
		result.push_back(1.0);
		for (size_t i = 0; i < number_of_nodes; i++) {
			if (i != node) {
				if (state.get_node_value(i) == 1) {
					result.push_back(1.0);
				}
				else {
					result.push_back(0.0);
				}
			}
		}
		return result;
	}
public:
	LikelihoodTester<Real_t>(size_t number, Real_t t):
		number_of_nodes{number},
		t_max{t} {}

	Real_t calculate_likelihood(const std::vector<Real_t> &beta, const std::vector<std::pair<State, Real_t>> &skeleton, const size_t node, const size_t node_past_value) {
		Real_t result = 0;
		for (size_t i = 0; i < skeleton.size() - 1; i++) {
			if (skeleton[i].first.get_node_value(node) == node_past_value) {
				auto intensity = get_intensity(beta, skeleton[i].first, node);
				if (skeleton[i + 1].first.get_node_value(node) != skeleton[i].first.get_node_value(node)) {
					result -= intensity;
				}
				result += std::exp(intensity + std::log(skeleton[i + 1].second - skeleton[i].second));
			}
		}
		if (skeleton.back().first.get_node_value(node) == node_past_value) {
			auto intensity = get_intensity(beta, skeleton.back().first, node);
			result += std::exp(intensity) * (t_max - skeleton.back().second);
		}
		return result / t_max;
	}

	std::vector<Real_t> calculate_gradient(const std::vector<Real_t> &beta, const std::vector<std::pair<State, Real_t>> &skeleton, const size_t node, const size_t node_past_value) {
		std::vector<Real_t> result;
		result.resize(beta.size());
		for (size_t i = 0; i < beta.size(); i++) {
			result[i] = 0.0;
		}
		for (size_t i = 0; i < skeleton.size() - 1; i++) {
			if (skeleton[i].first.get_node_value(node) == node_past_value) {
				auto intensity = get_intensity(beta, skeleton[i].first, node);
				auto predictive = get_predictive_vector(skeleton[i].first, node);
				if (skeleton[i + 1].first.get_node_value(node) != skeleton[i].first.get_node_value(node)) {
					for (size_t c = 0; c < result.size(); c++) {
						result[c] -= predictive[c];
					}
				}
				for (size_t c = 0; c < result.size(); c++) {
					result[c] += std::exp(intensity + std::log(skeleton[i + 1].second - skeleton[i].second)) * predictive[c];
				}
			}
		}
		if (skeleton.back().first.get_node_value(node) == node_past_value) {
			auto intensity = get_intensity(beta, skeleton.back().first, node);
			auto predictive = get_predictive_vector(skeleton.back().first, node);
			for (size_t c = 0; c < result.size(); c++) {
				result[c] += std::exp(intensity) * (t_max - skeleton.back().second) * predictive[c];
			}
		}
		for (size_t c = 0; c < result.size(); c++) {
			result[c] = result[c] / t_max;
		}
		return result;
	}

};
#endif // !LIKELIHOOD_TESTER_H
