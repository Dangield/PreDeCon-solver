#include "predecon.hpp"

int main(int argc, char const *argv[])
{
	predecon solver;
	std::string filename = "../test_data";
	solver.loadDataFromFile(filename);
	solver.printData();
	return 0;
}