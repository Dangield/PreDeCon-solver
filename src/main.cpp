#include "predecon.hpp"

int main(int argc, char const *argv[])
{
	predecon solver;
	std::string filename = "../test_data";
	solver.loadDataFromFile(filename);
	solver.setParameters(2.0, 0.0, 0, 0);
	solver.calculateENeighbours();
	// solver.printENeighbours();
	solver.printData();
	return 0;
}