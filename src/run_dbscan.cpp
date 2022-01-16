#include "dbscan.hpp"
#include <cstring>

int main(int argc, char const *argv[])
{
	dbscan solver;
	if (argc > 2) {
		for (int i = 2; i < argc; i++) {
			if (!strcmp(argv[i], "-e")) {
				solver.setEpsilon(std::stof(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-m")) {
				solver.setMinPts(std::stoi(argv[i+1]));
				i++;
			}
			if (!strcmp(argv[i], "-m_order")) {
				solver.setMiknowskiOrder(DistanceMetric(std::stoi(argv[i+1])));
				i++;
			}

		}
	}
	solver.run(argv[1]);
	return 0;
}