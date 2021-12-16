#include <string>
#include <vector>
#include "sample.hpp"

class predecon
{
private:
	float epsilon; // minimum distance threshold for e-neighbourhood
	float delta; // upper bound for the variance in an attribute
	int lambda; // maximum number of attributes that have a low variance
	int mi; // minimum threshold for weighted e-neighbourhood size
	float kappa; // constant for subspace preference vector
	std::vector<sample> data; // data for claustering
	std::vector<std::string> attribute_names; // names of attributes
public:
	predecon();
	~predecon();
	void setParameters(int, int, int, int, int);
	void printParameters();
	void loadDataFromFile(std::string);
	void printData();
};