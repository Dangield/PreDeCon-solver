#include "predecon.hpp"

int main(int argc, char const *argv[])
{
	predecon solver;
	std::string filename = "../test_data";
	solver.loadDataFromFile(filename);
	solver.setParameters(2, 0.3, 2, 2, 10.0);
	solver.calculateENeighbours();
	solver.printData();
	solver.calculateVariances();
	solver.printData();
	solver.calculateSubspacePreferenceVectors();
	solver.printData();
	solver.calculatePreferenceWeightedENeighbours();
	solver.printData();
	solver.solve();
	solver.printData();
	return 0;
}