#include "predecon.hpp"

int main(int argc, char const *argv[])
{
	predecon solver;
	std::string filename = "../test_data";
	solver.loadDataFromFile(filename);
	// solver.setParameters(2, 0.625, 2, 2, 10.0);
	solver.setParameters(3, 0.625, 2, 2, 10.0);
	// solver.setParameters(2, 0.3, 2, 2, 10.0);
	solver.setDistanceMetric(Minkowsky2);
	solver.calculateENeighbours();
	solver.printData(true);
	solver.calculateVariances();
	solver.printData(true);
	solver.calculateSubspacePreferenceVectors();
	solver.printData(true);
	solver.calculatePreferenceWeightedENeighbours();
	solver.solve();
	solver.printData();
	solver.plotData();
	return 0;
}