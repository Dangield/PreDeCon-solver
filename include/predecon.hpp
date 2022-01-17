#include <string>
#include <vector>
#include "sample.hpp"

enum DistanceMetric {MinkowskiInf, Minkowski1, Minkowski2};
// enum Property {Unknown, Core, Noise};
// static std::string property_strings[] = {"?", "Core", "Noise"};

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
	std::string filename;
	std::string dataset_name;
	int attribute_amount, sample_amount, current_cluster_id;
	int64_t start_time, file_read_time, calculate_e_neighbours_time,
		calculate_variances_time, calculate_subspace_preference_vectors_time,
		calculate_preference_weighted_e_neighbours_time, clustering_time,
		calculate_rand_time, calculate_purity_time, calculate_silhouette_coefficient_time,
		write_out_file_time, write_debug_file_time;
	DistanceMetric metric;
	std::vector<std::vector<int>> e_neighbours;
	std::vector<std::vector<float>> variances;
	std::vector<std::vector<float>> subspace_preference_vectors;
	std::vector<std::vector<int>> preference_weighted_e_neighbours;
	std::vector<int> num_of_distance_calc, point_type, cluster_id;
	std::vector<std::string> original_cluster_id;
	std::vector<int> original_int_cluster_id;
	int tp, tn, pairs;
	float rand, purity, silhouette_coefficient;
	// std::vector<Property> core;
	// std::vector<int> cluster;
	// std::vector<int> queue_theta, queue_R;
	bool use_TI;
	std::vector<float> reference_point;
	std::vector<float> reference_point_distance;
	std::vector<sample> sorted_data;
	int64_t calculate_reference_point_distance_time, sort_data_time;

	void readFile(std::string);
	void calculateENeighbours();
	void calculateENeighbours(int);
	float calculateDistance(int, int);
	void calculateVariances();
	void calculateSubspacePreferenceVectors();
	void calculatePreferenceWeightedENeighbours();
	void calculatePreferenceWeightedENeighbours(int);
	float calculateWeightedDistance(int, int);
	void performClustering();
	// float calculateDistance(sample, sample);
	// float calculateDistanceEuclidean(sample, sample);
	// float calculateDistanceMinkowsky1(sample, sample);
	// float calculateDistanceMinkowsky2(sample, sample);
	// float calculateDistanceMinkowskyInf(sample, sample);
	// std::vector<float> calculateVariances(sample, std::vector<std::string>);
	// std::vector<float> calculateSubspacePreferenceVectors(std::vector<float>);
	// std::vector<std::string> calculatePreferenceWeightedENeighbours(sample, std::vector<float>);
	// float calculateWeightedDistance(sample, sample, std::vector<float>);
	// float calculateWeightedDistanceEuclidean(sample, sample, std::vector<float>);
	// float calculateWeightedDistanceMinkowsky1(sample, sample, std::vector<float>);
	// float calculateWeightedDistanceMinkowsky2(sample, sample, std::vector<float>);
	// float calculateWeightedDistanceMinkowskyInf(sample, sample, std::vector<float>);
	bool isCore(int);
	// bool isNoise(int);
	int calculateSubspacePreferenceDimensionality(int);
	void writeOUTFile();
	void writeSTATFile();
	void writeDEBUGFile();
	void calculateRand();
	void calculatePurity();
	void calculateSilhouetteCoefficient();
	void calculateReferencePointDistances();
	float calculateDistanceToReferencePoint(int);
	void sortData();
	void calculateENeighboursTI();
public:
	predecon();
	~predecon();
	void setParameters(float, float, int, int, float);
	void setEpsilon(float);
	void setDelta(float);
	void setLambda(int);
	void setMi(int);
	void setKappa(float);
	void setMinkowskiOrder(DistanceMetric);
	void useTI(bool);
	// void printParameters();
	// void loadDataFromFile(std::string, bool dataHasIds = true);
	// void printData(bool waitForInput = false);
	// void setDistanceMetric(DistanceMetric);
	// void calculateENeighbours();
	// void calculateVariances();
	// void calculateSubspacePreferenceVectors();
	// void calculatePreferenceWeightedENeighbours();
	// void solve();
	// void plotData();
	void run(std::string);
};