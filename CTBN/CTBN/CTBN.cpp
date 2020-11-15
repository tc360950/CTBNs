// CTBN.cpp: definiuje punkt wejścia dla aplikacji.
//
#include <random>
#include <chrono>


#include "CTBN.h"
#include "bob_dylan.h"
#include "models/list_model.h"
using namespace std;
std::mt19937 generator;
bool random_bit() {
	std::uniform_int_distribution<int> distrib(0, 1);
	auto bit = distrib(generator);
	return  bit == 1 ? true : false;
}
int main()
{	
	std::vector<bool> b;
	std::vector<double> dupa;
	std::vector<double> dupa2;
	for (int i = 0; i < 20; i++) {
		b.push_back(random_bit());
		dupa2.push_back(b[i]);
		dupa.push_back(i);
	}

	{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		double r1 = 0.0;
		for (int i = 0; i < 500000; i++) {
			for (int i = 0; i < 20; i++) {
					r1 += b[i] * dupa[i];
			}
		}
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
	}
	{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		double r2 = 0.0;
		for (int i = 0; i < 500000; i++) {
			for (int i = 0; i < 20; i++) {
				r2 += dupa[i] * dupa2[i];
			}
		}
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
	}

	cout << "Witaj, CMake." << endl;
	BobDylan<double, ListModel<double>> bob;
	bob.simulate(20, 12314214, 50);
	return 0;
}
