#include <string>
#include <vector>
#include "sample.hpp"

enum DistanceMetric {Euclidean, Minkowsky1, Minkowsky2, MinkowskyInf};
enum Property {Unknown, Core, Noise};
static std::string property_strings[] = {"?", "Core", "Noise"};

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
	int attribute_amount;
	DistanceMetric distance_metric;
	std::vector<std::vector<std::string>> e_neighbours;
	std::vector<std::vector<float>> variances;
	std::vector<std::vector<float>> subspace_preference_vectors;
	std::vector<std::vector<std::string>> preference_weighted_e_neighbours;
	std::vector<Property> core;
	std::vector<int> cluster;
	int current_cluster_id;
	std::vector<int> queue_theta, queue_R;

	std::vector<std::string> calculateENeighbours(sample);
	float calculateDistance(sample, sample);
	float calculateDistanceEuclidean(sample, sample);
	float calculateDistanceMinkowsky1(sample, sample);
	float calculateDistanceMinkowsky2(sample, sample);
	float calculateDistanceMinkowskyInf(sample, sample);
	std::vector<float> calculateVariances(sample, std::vector<std::string>);
	std::vector<float> calculateSubspacePreferenceVectors(std::vector<float>);
	std::vector<std::string> calculatePreferenceWeightedENeighbours(sample, std::vector<float>);
	float calculateWeightedDistance(sample, sample, std::vector<float>);
	float calculateWeightedDistanceEuclidean(sample, sample, std::vector<float>);
	float calculateWeightedDistanceMinkowsky1(sample, sample, std::vector<float>);
	float calculateWeightedDistanceMinkowsky2(sample, sample, std::vector<float>);
	float calculateWeightedDistanceMinkowskyInf(sample, sample, std::vector<float>);
	bool isCore(int);
	bool isNoise(int);
	int calculateSubspacePreferenceDimensionality(int);
public:
	predecon();
	~predecon();
	void setParameters(float, float, int, int, float k = 100);
	void printParameters();
	void loadDataFromFile(std::string, bool dataHasIds = true);
	void printData(bool waitForInput = false);
	void setDistanceMetric(DistanceMetric);
	void calculateENeighbours();
	void calculateVariances();
	void calculateSubspacePreferenceVectors();
	void calculatePreferenceWeightedENeighbours();
	void solve();
	void plotData();
};