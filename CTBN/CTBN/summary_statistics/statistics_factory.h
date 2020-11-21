#ifndef STATISTICS_FACTORY_H
#define STATISTICS_FACTORY_H
#include <set>
#include <utility>

#include "statistics.h"
#include "../models/model_data.h"
#include "../bob_dylan.h"
#include "../utils/result.h"
//TODO trzbea to torche inaczje robic przy interakcjach 
template <class Real_t, class Model> class StatisticsFactory {

	std::set<std::pair<size_t, size_t>> gather_edges(const size_t number_of_nodes, SimulationResult<Real_t> simulation_result) {
		std::set<std::pair<size_t, size_t>> result;
		for (size_t node = 0; node < number_of_nodes; node++) {
			for (size_t i = 1; i < simulation_result.inference_results[2 * node].beta.size(); i++) {
				if (simulation_result.inference_results[2 * node].beta[i] != 0.0 || simulation_result.inference_results[2 * node + 1].beta[i] != 0.0) {
					size_t node2 = i - 1 >= node ? i : i - 1;
					result.insert(std::make_pair(node2, node));
				}
			}
		}
		return result;
	}

	std::set<std::pair<size_t, size_t>> gather_real_edges(const size_t number_of_nodes, ModelData<Real_t> model_data) {
		std::set<std::pair<size_t, size_t>> result;
		for (size_t i = 0; i < number_of_nodes; i++) {
			for (size_t j = 0; j < number_of_nodes; j++) {
				if (model_data.dependence_structure[i][j]) {
					result.insert(std::make_pair(i, j));
				}
			}
		}
		return result;
	}

	Real_t calculate_power(const std::set<std::pair<size_t, size_t>> &inferred_edges, const std::set<std::pair<size_t, size_t>> &real_edges) {
		Real_t selected = 0.0;
		for (auto x : inferred_edges) {
			if (real_edges.find(x) != real_edges.end()) {
				selected++;
			}
		}
		return selected / (Real_t) real_edges.size();
	}

	Real_t calculate_FDR(const std::set<std::pair<size_t, size_t>> &inferred_edges, const std::set<std::pair<size_t, size_t>> &real_edges) {
		Real_t selected = 0.0;
		for (auto x : inferred_edges) {
			if (real_edges.find(x) == real_edges.end()) {
				selected++;
			}
		}
		return selected / (Real_t)inferred_edges.size();
	}
public:
	Statistics<Real_t> convert(const size_t number_of_nodes, const Real_t time, ModelData<Real_t> model_data, SimulationResult<Real_t> simulation_result) {
		auto inferred_edges = gather_edges(number_of_nodes, simulation_result);
        for (auto edge : inferred_edges) {
            std::cout << edge.first << " " << edge.second << "\n";

        }
		auto real_edges = gather_real_edges(number_of_nodes, model_data);
		auto power = calculate_power(inferred_edges, real_edges);
		auto FDR = calculate_FDR(inferred_edges, real_edges);
		Real_t MD = (Real_t)inferred_edges.size();
		return Statistics<Real_t>{number_of_nodes, time, power, FDR, MD};
	}
};


#endif // !STATISTICS_FACTORY_H
