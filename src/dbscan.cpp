#include "dbscan.hpp"
#include <fstream>
#include <iostream>
#include <cstring>
#include <chrono>
#include <ctime>
#include <sys/time.h>
#include <cmath>
#include <algorithm>

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::system_clock;

dbscan::dbscan() {
	epsilon = 0;
	minPts = 0;
	metric = Minkowski2;
	printf("DBSCAN initialized\n");
}

dbscan::~dbscan() {

}

void dbscan::setParameters(float e, int m) {
	epsilon = e;
	minPts = m;
	printf("Setting DBSCAN parameters: eps %f, minPts %d\n", epsilon, minPts);
}

void dbscan::setEpsilon(float e) {
	epsilon = e;
	printf("Setting DBSCAN parameters: eps %f\n", epsilon);
}

void dbscan::setMinPts(int m) {
	minPts = m;
	printf("Setting DBSCAN parameters: minPts %d\n", minPts);
}

void dbscan::setMinkowskiOrder(DistanceMetric m) {
	metric = m;
	if (metric == 1) printf("Setting distance metric to: Minkowski1\n");
	else if (metric == 2) printf("Setting distance metric to: Minkowski2\n");
	else if (metric == 0) printf("Setting distance metric to: MinkowskiInf\n");
	else printf("Warning: unknown distance metric set\n");
}

void dbscan::readFile(std::string fname){
	filename = fname;
	std::ifstream file(filename);
	if (!file) {
		printf("File with data doesn't exist.\n");
		return;
	}
	std::string line, value;
	size_t pos = 0;
	bool data_with_class_names = false;

	while(getline(file, line)) {
		if (line[0] == '%') continue;
		if (line[0] == '@') {
			pos = line.find(std::string(" "));
			value = line.substr(1, pos-1);
			line.erase(0, pos + 1);
			if (!strcmp(value.c_str(), "relation") or !strcmp(value.c_str(), "RELATION")){
				dataset_name = line;
			} else if (!strcmp(value.c_str(), "attribute") or !strcmp(value.c_str(), "ATTRIBUTE")){
				pos = line.find(std::string(" "));
				value = line.substr(0, pos);
				line.erase(0, pos + 1);
				if (!strcmp(value.c_str(), "class") || !strcmp(value.c_str(), "CLASS")) {
					data_with_class_names = true;
					line.erase(0, 1);
					while ((pos = line.find(',')) != std::string::npos) {
						pos = line.find(std::string(","));
						value = line.substr(0, pos);
						while(value[0] == ' ') value.erase(0, 1);
						if (strcmp(value.c_str(), "noise")) original_cluster_id.push_back(value.substr(0, pos));
						line.erase(0, pos + 1);
					}
					pos = line.find(std::string("}"));
					value = line.substr(0, pos);
					while(value[0] == ' ') value.erase(0, 1);
					if (strcmp(value.c_str(), "noise")) original_cluster_id.push_back(value.substr(0, pos));
					continue;
				}
				else {
					attribute_names.push_back(std::string(value));
				}
			} else if (!strcmp(value.c_str(), "data") || !strcmp(value.c_str(), "DATA")) {
				break;
			}
		}
	}
	// printf("Num of original classes: %lu\n", original_cluster_id.size());
	// for (int i = 0; i < original_cluster_id.size(); i++) printf("%s", original_cluster_id[i].c_str());

	attribute_amount = attribute_names.size();
	while(getline(file, line)) {
		if (line[0] == '%') continue;
		if (std::isdigit(line[0])) {
			sample s;
			while ((pos = line.find(',')) != std::string::npos) {
				value = line.substr(0, pos);
				s.pushAttribute(std::stof(value));
				line.erase(0, pos + 1);
			}
			if (!data_with_class_names) s.pushAttribute(std::stof(line));
			else {
				if (!strcmp(line.c_str(), "noise")) original_int_cluster_id.push_back(-1);
				else {
					for (int i = 0; i < original_cluster_id.size(); i++) {
						// printf("%i", i);
						if (!strcmp(line.c_str(), original_cluster_id[i].c_str())) {
							original_int_cluster_id.push_back(i + 1);
							// printf("\n");
							break;
						}
					}
				}
			}
			data.push_back(s);
			// printf("%i", original_int_cluster_id[0]);
		}
	}
	sample_amount = data.size();
	// for (int i = 0; i < sample_amount; i++) std::cout << original_int_cluster_id[i] << std::endl;
	// for (int i = 0; i < sample_amount; i++) printf("%d\t%d\n", data[i].getId(), original_int_cluster_id[i]);
	// printf("%s\n", dataset_name.c_str());
	// for (std::string a : attribute_names) printf("%s\n", a.c_str());
	// for (sample s : data) for (float f : s.getAttributes()) printf("%f\n", f);
	printf("Dataset loaded: %s\n", dataset_name.c_str());
}

void dbscan::run(std::string fname) {
	start_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	// printf("%ld\n",start_time);
	readFile(fname);
	file_read_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",file_read_time - start_time);
	calculateENeighbours();
	calculate_e_neighbours_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_e_neighbours_time - file_read_time);
	performClustering();
	clustering_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",clustering_time - calculate_e_neighbours_time);
	// for (int i = 0; i < sample_amount; i++) printf("%i\t%i\n", cluster_id[i], point_type[i]);
	// remapClusterIds();
	calculateRand();
	calculate_rand_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_rand_time - clustering_time);
	calculatePurity();
	calculate_purity_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_purity_time - calculate_rand_time);
	calculateSilhouetteCoefficient();
	calculate_silhouette_coefficient_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_silhouette_coefficient_time - calculate_purity_time);
	writeOUTFile();
	write_out_file_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",write_out_file_time - calculate_silhouette_coefficient_time);
	writeDEBUGFile();
	write_debug_file_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",write_debug_file_time - write_out_file_time);
	printf("Total runtime: %ld\n", write_debug_file_time - start_time);
	writeSTATFile();
}

void dbscan::calculateENeighbours() {
	for (int i = 0; i < sample_amount; i++) {
		e_neighbours.push_back(std::vector<int>(1, i));
		num_of_distance_calc.push_back(0);
	}
	for (int i = 0; i < sample_amount; i++) {
		std::cout << i << '\r';
		calculateENeighbours(i);
	}
	printf("Calculated neighbours for the data\n");
	// for (std::vector<int> v : e_neighbours) {
	// 	for (int i : v) printf("%d ", i);
	// 	printf("\n");
	// }
	// for (int i : num_of_distance_calc) printf("%d\n", i);
}

void dbscan::calculateENeighbours(int p) {
	for (int i = p + 1; i < sample_amount; i++) {
		num_of_distance_calc[p]++;
		num_of_distance_calc[i]++;
		if (calculateDistance(p, i) <= epsilon) {
			e_neighbours[p].push_back(i);
			e_neighbours[i].push_back(p);
		}
	}
}

float dbscan::calculateDistance(int p, int q) {
	float distance = 0;
	if (metric == Minkowski2){
		for (int i = 0; i < attribute_amount; i++) distance += pow(data[p].getAttribute(i) - data[q].getAttribute(i), 2.0);
		return pow(distance, 0.5);
	} else if (metric == Minkowski1) {
		for (int i = 0; i < attribute_amount; i++) distance += abs(data[p].getAttribute(i) - data[q].getAttribute(i));
		return distance;
	} else if (metric == MinkowskiInf) {
		for (int i = 0; i < attribute_amount; i++) distance = std::max(std::abs(data[p].getAttribute(i) - data[q].getAttribute(i)), distance);
		return distance;
	}
}

bool dbscan::isCore(int i) {
	return e_neighbours[i].size() >= minPts;
}

void dbscan::performClustering() {
	current_cluster_id = 0;
	for (int i = 0; i < sample_amount; i++) {
		point_type.push_back(0);
		cluster_id.push_back(0);
	}
	std::vector<int> seeds;
	for (int i = 0; i < sample_amount; i++) {
		if (cluster_id[i] == 0) {
			if (isCore(i)) {
				point_type[i] = 1;
				current_cluster_id++;
				cluster_id[i] = current_cluster_id;
				for (int n : e_neighbours[i]) {
					if (cluster_id[n] == 0) {
						cluster_id[n] = current_cluster_id;
						seeds.push_back(n);
					} else if (point_type[n] == -1) {
						cluster_id[n] = current_cluster_id;
						point_type[n] = 0;
					}
				}
				// for (int i = 0; i < sample_amount; i++) printf("%i\t%i\n", cluster_id[i], point_type[i]);
				// std::cin.get();
				while (!seeds.empty()) {
					if (isCore(seeds[0])) {
						point_type[seeds[0]] = 1;
						for (int n : e_neighbours[seeds[0]]) {
							if (cluster_id[n] == 0) {
								cluster_id[n] = current_cluster_id;
								seeds.push_back(n);
							} else if (point_type[n] == -1) {
								cluster_id[n] = current_cluster_id;
								point_type[n] = 0;
							}
						}
					}
					seeds.erase(seeds.begin());
				}
			} else {
				point_type[i] = -1;
				cluster_id[i] = -1;
			}
			// for (int i = 0; i < sample_amount; i++) printf("%i\t%i\n", cluster_id[i], point_type[i]);
			// std::cin.get();
		}
	}
	printf("Clustering complete\n");
}

void dbscan::writeOUTFile() {
	std::ofstream file;
	std::string out_fname("OUT_DBSCAN_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_m" + std::to_string(minPts) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + ".csv");
	file.open(out_fname);
	if (file.is_open()){
		file << "point id, ";
		for (std::string att : attribute_names) file << att + ", ";
		file << "# of distance/similarity calculations, point type, CId\n";
		for (int i = 0; i < sample_amount; i++){
			file << data[i].getId() << ", ";
			for (float att : data[i].getAttributes()) file << att << ", ";
			file << num_of_distance_calc[i] << ", " << point_type[i] << ", " << cluster_id[i] << "\n";
		}
		file.close();
	}
	printf("OUT file created\n");
}

void dbscan::writeDEBUGFile() {
	std::ofstream file;
	std::string out_fname("DEBUG_DBSCAN_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_m" + std::to_string(minPts) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + ".csv");
	file.open(out_fname);
	if (file.is_open()){
		file << "point id, size of eps_neighbourhood\n";
		for (int i = 0; i < sample_amount; i++){
			file << data[i].getId() << ", " << e_neighbours[i].size() << "\n";
		}
		file.close();
	}
	printf("DEBUG file created\n");
}

void dbscan::writeSTATFile() {
	std::ofstream file;
	std::string out_fname("STAT_DBSCAN_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_m" + std::to_string(minPts) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + ".csv");
	file.open(out_fname);
	if (file.is_open()){
		file << "Input file: " << filename << "\n";
		file << "Number of dimensions of a point: " << attribute_amount << "\n";
		file << "Number of points: " << sample_amount << "\n";
		file << "Reading a file runtime[us]: " << file_read_time - start_time << "\n";
		file << "Calculating epsilon neighbours runtime[us]: " << calculate_e_neighbours_time - file_read_time << "\n";
		file << "Clustering runtime[us]: " << clustering_time - calculate_e_neighbours_time << "\n";
		file << "Calculating RAND measure runtime[us]: " << calculate_rand_time - clustering_time << "\n";
		file << "Calculating Purity measure runtime[us]: " << calculate_purity_time - calculate_rand_time << "\n";
		file << "Calculating Silhouette Coefficient runtime[us]: " << calculate_silhouette_coefficient_time - calculate_purity_time << "\n";
		file << "Writing to OUT file runtime[us]: " << write_out_file_time - calculate_silhouette_coefficient_time << "\n";
		file << "Writing to DEBUG file runtime[us]: " << write_debug_file_time - write_out_file_time << "\n";
		file << "Total runtime [us]: " << write_debug_file_time - start_time << "\n";
		file << "Number of discovered clusters: " << current_cluster_id << "\n";
		int n = 0, c = 0, b = 0;
		float avg = 0;
		for (int i = 0; i < sample_amount; i++) {
			if (point_type[i] == -1) n++;
			else if (point_type[i] == 0) b++;
			else if (point_type[i] == 1) c++;
			avg += num_of_distance_calc[i];
		}
		avg /= sample_amount;
		file << "Number of discovered noise points: " << n << "\n";
		file << "Number of discovered core points: " << c << "\n";
		file << "Number of discovered border points: " << b << "\n";
		file << "Average number of distance calculations of a point to other points in the data set: " << avg << "\n";
		file << "TP: " << tp << "\n";
		file << "TN: " << tn << "\n";
		file << "NUmber of pairs of points: " << pairs << "\n";
		file << "RAND: " << rand << "\n";
		file << "Purity: " << purity << "\n";
		file << "Silhouette Coefficient: " << silhouette_coefficient << "\n";
		file.close();
	}
	printf("STAT file created\n");
}

void dbscan::remapClusterIds() {
	std::vector<int> samples_per_original_class(original_cluster_id.size(), 0);
	for (int i = 0; i < sample_amount; i++) {
		if (original_int_cluster_id[i] == -1) continue;
		else samples_per_original_class[original_int_cluster_id[i] - 1]++;
	}
	// for (int i : samples_per_original_class) printf("%d", i);
	std::vector<int> map_to(current_cluster_id + 1, 0);
	while (samples_per_original_class[std::distance(samples_per_original_class.begin(), std::max_element(samples_per_original_class.begin(), samples_per_original_class.end()))] > 0){
		int max_pos = std::distance(samples_per_original_class.begin(), std::max_element(samples_per_original_class.begin(), samples_per_original_class.end()));
		// printf("Original cid: %d with sample amount: %d\n", max_pos + 1, samples_per_original_class[max_pos]);
		samples_per_original_class[max_pos] = 0;
		std::vector<int> calc_clusters_for_this_cluster(current_cluster_id, 0);
		for (int i = 0; i < sample_amount; i++) {
			if (original_int_cluster_id[i] == max_pos + 1) {
				if (cluster_id[i] != -1) calc_clusters_for_this_cluster[cluster_id[i] - 1]++;
			}
		}
		map_to[std::distance(calc_clusters_for_this_cluster.begin(), std::max_element(calc_clusters_for_this_cluster.begin(), calc_clusters_for_this_cluster.end())) + 1] = max_pos + 1;
		// for (int i : calc_clusters_for_this_cluster) printf("%d ", i);
		// printf("\n");
	}
	// for (int i : map_to) printf("%d ", i);
	// printf("\n");
	if (current_cluster_id > original_cluster_id.size()) {
		int next_id = original_cluster_id.size() + 1;
		for (int i = 1; i < current_cluster_id + 1; i++) {
			if (map_to[i] == 0) {
				map_to[i] = next_id;
				next_id++;
			}
		}
	}
	// for (int i : map_to) printf("%d ", i);
	// printf("\n");
	for (int i = 0; i < sample_amount; i++) cluster_id[i] = map_to[cluster_id[i]];
}

void dbscan::calculateRand() {
	pairs = sample_amount*(sample_amount-1)/2;
	// calculate TP and TN
	tp = 0;
	tn = 0;
	for (int p = 0; p < sample_amount; p++) {
		for (int q = p; q < sample_amount; q++) {
			if (original_int_cluster_id[p] == original_int_cluster_id[q] && cluster_id[p] == cluster_id[q]) tp++;
			if (original_int_cluster_id[p] != original_int_cluster_id[q] && cluster_id[p] != cluster_id[q]) tn++;
		}
	}
	rand = (tp + tn)/pairs;
	printf("Calculated RAND\n");
}

void dbscan::calculatePurity(){
	purity = 0;
	float one_cluster_purity, partial_purity;
	std::vector<int> samples_per_class(original_cluster_id.size()+1, 0);
	for (int i = 0; i < current_cluster_id + 1; i++){
		one_cluster_purity = 0;
		for (int j = 0; j < original_cluster_id.size() + 1; j++) {
			partial_purity = 0;
			for (int s = 0; s < sample_amount; s++) {
				if (i == 0 && j == 0 && cluster_id[s] == -1 && original_int_cluster_id[s] == -1) partial_purity++;
				else if (i == 0 && cluster_id[s] == -1 && original_int_cluster_id[s] == j) partial_purity++;
				else if (j == 0 && cluster_id[s] == i && original_int_cluster_id[s] == -1) partial_purity++;
				else if (cluster_id[s] == i && original_int_cluster_id[s] == j) partial_purity++;
			}
			one_cluster_purity = std::max(one_cluster_purity, partial_purity);
		}
		purity += one_cluster_purity/sample_amount;
	}
	printf("Calculated Purity\n");
}

void dbscan::calculateSilhouetteCoefficient() {
	silhouette_coefficient = 0;
	std::vector<int> num_of_points_in_cluster(current_cluster_id + 1, 0);
	for (int i = 0; i < sample_amount; i++) {
		if (cluster_id[i] != -1) num_of_points_in_cluster[cluster_id[i]]++;
	}
	float a, b, b_cluster, s;
	for (int i = 0; i < sample_amount; i++) {
		if (cluster_id[i] == -1) continue;
		a = 0;
		b = 0;
		// calculate a
		for (int j = 0; j < sample_amount; j++) {
			if (i == j) continue;
			if (cluster_id[i] == cluster_id[j]) a += calculateDistance(i, j)/(num_of_points_in_cluster[cluster_id[i]] - 1);
		}
		// calculate b
		for (int c = 1; c < original_cluster_id.size() + 1; c++) {
			if (c == cluster_id[i]) continue;
			b_cluster = 0;
			for (int j = 0; j < sample_amount; j++) {
				if (cluster_id[j] == c) b_cluster += calculateDistance(i, j)/num_of_points_in_cluster[c];
			}
			if (b == 0) b = b_cluster;
			else b = std::min(b, b_cluster);
		}
		s = (b - a)/std::max(a, b);
		silhouette_coefficient = std::max(silhouette_coefficient, s);
	}
	printf("Calculated Silhouette Coefficient\n");
}