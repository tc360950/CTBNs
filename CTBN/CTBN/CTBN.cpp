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
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	cout << "Witaj, CMake." << endl;
	BobDylan<double, ListModel<double>> bob;
	bob.simulate(20, 12314214, 50);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[µs]" << std::endl;
	
	while (true) {}
	return 0;
}
