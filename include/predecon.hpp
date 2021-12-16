#include <string>
#include <vector>
#include "sample.hpp"

enum DistanceMetric {Euclidean, Minkowsky1, Minkowsky2, MinkowskyInf};

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
	DistanceMetric distance_metric;
	int attribute_amount;
	std::vector<std::vector<std::string>> e_neighbours;
	float calculateDistanceEuclidean(sample, sample);
	std::vector<std::string> calculateENeighbours(sample);
public:
	predecon();
	~predecon();
	void setParameters(float, float, int, int, float k = 100);
	void printParameters();
	void loadDataFromFile(std::string);
	void printData();
	float calculateDistance(sample, sample);
	void setDistanceMetric(DistanceMetric);
	void calculateENeighbours();
	void printENeighbours();
};