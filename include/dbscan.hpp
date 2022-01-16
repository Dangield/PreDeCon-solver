#include <string>
#include <vector>
#include "sample.hpp"

enum DistanceMetric {MinkowskiInf, Minkowski1, Minkowski2};

class dbscan {
	private:
		float epsilon;
		int minPts;
		DistanceMetric metric;
		std::vector<sample> data; // data for claustering
		std::vector<std::string> attribute_names; // names of attributes
		std::string filename;
		std::string dataset_name;
		int attribute_amount, sample_amount, current_cluster_id;
		int64_t start_time, file_read_time, calculate_e_neighbours_time, clustering_time, calculate_rand_time, calculate_purity_time, calculate_silhouette_coefficient_time, write_out_file_time, write_debug_file_time;
		std::vector<std::vector<int>> e_neighbours;
		std::vector<int> num_of_distance_calc, point_type, cluster_id;
		std::vector<std::string> original_cluster_id;
		std::vector<int> original_int_cluster_id;
		int tp, tn, pairs;
		float rand, purity, silhouette_coefficient;

		void readFile(std::string);
		void calculateENeighbours();
		void calculateENeighbours(int);
		float calculateDistance(int, int);
		void performClustering();
		bool isCore(int);
		void writeOUTFile();
		void writeSTATFile();
		void writeDEBUGFile();
		void remapClusterIds();
		void calculateRand();
		void calculatePurity();
		void calculateSilhouetteCoefficient();
	public:
		dbscan();
		~dbscan();
		void setParameters(float, int);
		void setEpsilon(float);
		void setMinPts(int);
		void setMiknowskiOrder(DistanceMetric);
		void run(std::string);
};