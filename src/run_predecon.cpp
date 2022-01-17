#include "predecon.hpp"
#include <cstring>

int main(int argc, char const *argv[])
{
	predecon solver;
	if (argc > 2) {
		for (int i = 2; i < argc; i++) {
			if (!strcmp(argv[i], "-e")) {
				solver.setEpsilon(std::stof(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-d")) {
				solver.setDelta(std::stof(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-l")) {
				solver.setLambda(std::stoi(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-m")) {
				solver.setMi(std::stoi(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-k")) {
				solver.setKappa(std::stof(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-m_order")) {
				solver.setMinkowskiOrder(DistanceMetric(std::stoi(argv[i+1])));
				i++;
			}

		}
	}
	solver.run(argv[1]);
	return 0;
	// predecon solver;
	// std::string filename = "../test_data";
	// solver.loadDataFromFile(filename);
	// // solver.setParameters(2, 0.625, 2, 2, 10.0);
	// solver.setParameters(3, 0.625, 2, 2, 10.0);
	// solver.setParameters(3, 0.625, 2, 4, 10.0);
	// // solver.setParameters(2, 0.3, 2, 2, 10.0);
	// solver.setDistanceMetric(Minkowsky2);
	// solver.calculateENeighbours();
	// solver.printData(true);
	// solver.calculateVariances();
	// solver.printData(true);
	// solver.calculateSubspacePreferenceVectors();
	// solver.printData(true);
	// solver.calculatePreferenceWeightedENeighbours();
	// solver.solve();
	// solver.printData();
	// solver.plotData();
	// return 0;
}