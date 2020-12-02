#ifndef RANDOM_TREE_H
#define RANDOM_TREE_H

#include <random>
#include <vector>

template <class Real_t> class RandomTree {
private:
	std::mt19937 generator;


	void generate_random_tree(std::vector<std::vector<bool>> &result_holder) {
		std::vector<bool> used = result_holder[0];
		std::fill(used.begin(), used.end(), false);
		size_t used_count = 1;
		size_t current_node = random_size_t(used.size());
		used[current_node] = true;
		while (used_count < used.size()) {
			size_t node = random_size_t_except(used.size(), current_node);
			if (!used[node]) {
				result_holder[node][current_node] = true;
				used[node] = true;
				used_count++;
			}
			current_node = node;
		}
	}

	void generate_random_binary_tree(std::vector<std::vector<bool>> &result_holder) {
		std::vector<bool> used = result_holder[0];
		std::fill(used.begin(), used.end(), false);
		size_t used_count = 1;
		size_t current_node = random_size_t(used.size());
		used[current_node] = true;
		while (used_count < used.size()) {
			while (has_degree_2(result_holder, current_node) || !used[current_node]) {
				current_node = random_size_t_except(used.size(), current_node);
			}
			size_t node = random_size_t_except(used.size(), current_node);
			if (!used[node]) {
				result_holder[node][current_node] = true;
				used[node] = true;
				used_count++;
			}
			current_node = node;
		}
	}

	bool has_degree_2(std::vector<std::vector<bool>> &tree, size_t node) {
		size_t count = 0;
		for (size_t i = 0; i < tree[0].size(); i++) {
			if (tree[i][node]) {
				count++;
			}
		}
		return count == 2;
	}
	size_t random_size_t_except(size_t bound, size_t except) {
		size_t res = random_size_t(bound);
		while (res == except) {
			res = random_size_t(bound);
		}
		return res;
	}

	size_t random_size_t(size_t bound) {
		std::uniform_int_distribution<size_t> distrib(0, bound - 1);
		return distrib(generator);
	}

public:
	RandomTree<Real_t>(long seed) : generator{ seed } {

	}

	std::vector<std::vector<bool>> generate_random_tree(size_t nodes) {
		std::vector<std::vector<bool>> result;
		std::vector<bool> false_vec;
		false_vec.resize(nodes);
		std::fill(false_vec.begin(), false_vec.end(), false);
		for (size_t i = 0; i < nodes; i++) {
			result.push_back(false_vec);
		}
		generate_random_tree(result);
		return result;
	}

	std::vector<std::vector<bool>> generate_random_binary_tree(size_t nodes) {
		std::vector<std::vector<bool>> result;
		std::vector<bool> false_vec;
		false_vec.resize(nodes);
		std::fill(false_vec.begin(), false_vec.end(), false);
		for (size_t i = 0; i < nodes; i++) {
			result.push_back(false_vec);
		}
		generate_random_binary_tree(result);
		return result;
	}
};









#endif //RANDOM_TREE_H
