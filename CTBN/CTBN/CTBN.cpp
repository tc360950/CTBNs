// CTBN.cpp: definiuje punkt wejścia dla aplikacji.
//

#include "CTBN.h"
#include "bob_dylan.h"
#include "models/list_model.h"
using namespace std;

int main()
{
	cout << "Witaj, CMake." << endl;
	BobDylan<double, ListModel<double>> bob;
	bob.simulate(20, 12314214);
	return 0;
}
