#include "predecon.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <sys/time.h>
// #include "matplotlibcpp.h"

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::system_clock;

predecon::predecon() {
	epsilon = 0;
	delta = 0;
	lambda = 0;
	mi = 0;
	kappa = 0;
	metric = Minkowski2;
	use_TI = false;
	printf("Predecon solver initialized\n");
}

predecon::~predecon() {
}

void predecon::run(std::string fname) {
	start_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	readFile(fname);
	file_read_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",file_read_time - start_time);
	if (use_TI) {
		reference_point = std::vector<float>(attribute_amount, 0);
		calculateReferencePointDistances();
		calculate_reference_point_distance_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
		printf("Time spent: %ld\n",calculate_reference_point_distance_time - file_read_time);
		sortData();
		sort_data_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
		printf("Time spent: %ld\n",sort_data_time - calculate_reference_point_distance_time);
		calculateENeighboursTI();
		calculate_e_neighbours_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
		printf("Time spent: %ld\n",calculate_e_neighbours_time - sort_data_time);
	} else {
		calculateENeighbours();
		calculate_e_neighbours_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
		printf("Time spent: %ld\n",calculate_e_neighbours_time - file_read_time);
	}
	calculateVariances();
	calculate_variances_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_variances_time - calculate_e_neighbours_time);
	calculateSubspacePreferenceVectors();
	calculate_subspace_preference_vectors_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_subspace_preference_vectors_time - calculate_variances_time);
	calculatePreferenceWeightedENeighbours();
	calculate_preference_weighted_e_neighbours_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",calculate_preference_weighted_e_neighbours_time - calculate_subspace_preference_vectors_time);
	performClustering();
	clustering_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	printf("Time spent: %ld\n",clustering_time - calculate_e_neighbours_time);
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

void predecon::setParameters(float e, float d, int l, int m, float k) {
	epsilon = e;
	delta = d;
	lambda = l;
	mi = m;
	kappa = k;
	printf("Setting PREDECON parameters: epsilon %f, delta %f, lambda %i, mi %i, kappa %f\n", epsilon, delta, lambda, mi, kappa);
}

void predecon::setEpsilon(float e) {
	epsilon = e;
	printf("Setting PREDECON parameters: epsilon %f\n", epsilon);
}

void predecon::setDelta(float d) {
	delta = d;
	printf("Setting PREDECON parameters: delta %f\n", delta);
}

void predecon::setLambda(int l) {
	lambda = l;
	printf("Setting PREDECON parameters: lambda %d\n", lambda);
}

void predecon::setMi(int m) {
	mi = m;
	printf("Setting PREDECON parameters: mi %d\n", mi);
}

void predecon::setKappa(float k) {
	kappa = k;
	printf("Setting PREDECON parameters: kappa %f\n", kappa);
}

void predecon::setMinkowskiOrder(DistanceMetric m) {
	metric = m;
	if (metric == 1) printf("Setting distance metric to: Minkowski1\n");
	else if (metric == 2) printf("Setting distance metric to: Minkowski2\n");
	else if (metric == 0) printf("Setting distance metric to: MinkowskiInf\n");
	else printf("Warning: unknown distance metric set\n");
}

void predecon::useTI(bool ti) {
	use_TI = ti;
}

void predecon::readFile(std::string fname){
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


void predecon::calculateENeighbours() {
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

void predecon::calculateENeighbours(int p) {
	for (int i = p + 1; i < sample_amount; i++) {
		num_of_distance_calc[p]++;
		num_of_distance_calc[i]++;
		if (calculateDistance(p, i) <= epsilon) {
			e_neighbours[p].push_back(i);
			e_neighbours[i].push_back(p);
		}
	}
}

float predecon::calculateDistance(int p, int q) {
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

void predecon::calculateReferencePointDistances() {
	for (int i = 0; i < sample_amount; i++) {
		reference_point_distance.push_back(calculateDistanceToReferencePoint(i));
	}
	printf("Calculated distance to reference point\n");
}

float predecon::calculateDistanceToReferencePoint(int p) {
	float distance = 0;
	if (metric == Minkowski2){
		for (int i = 0; i < attribute_amount; i++) distance += pow(data[p].getAttribute(i) - reference_point[i], 2.0);
		return pow(distance, 0.5);
	} else if (metric == Minkowski1) {
		for (int i = 0; i < attribute_amount; i++) distance += abs(data[p].getAttribute(i) - reference_point[i]);
		return distance;
	} else if (metric == MinkowskiInf) {
		for (int i = 0; i < attribute_amount; i++) distance = std::max(std::abs(data[p].getAttribute(i) - reference_point[i]), distance);
		return distance;
	}
}

void predecon::sortData() {
	std::vector<float> temp = reference_point_distance;
	std::sort(temp.begin(), temp.end());
	for (float d : temp) {
		for (int i = 0; i < sample_amount; i++) {
			if (reference_point_distance[i] == d) {
				if (sorted_data.size() == 0) sorted_data.push_back(data[i]);
				else {
					bool s_in_ds = false;
					for (sample s : sorted_data) {
						if (s.getId() == i) {
							s_in_ds = true;
							break;
						}
					}
					if (s_in_ds) continue;
					else (sorted_data.push_back(data[i]));
				}
			}
		}
	}
	// for (int i = 0; i < sample_amount; i++) printf("%d %f\n", sorted_data[i].getId(), reference_point_distance[sorted_data[i].getId()]);
	printf("Data sorted due to reference point distance\n");
}

void predecon::calculateENeighboursTI() {
	for (int i = 0; i < sample_amount; i++) {
		e_neighbours.push_back(std::vector<int>(1, i));
		num_of_distance_calc.push_back(0);
	}
	for (int p = 0; p < sample_amount; p++) {
		std::cout << p << '\r';
		for (int q = p + 1; q < sample_amount; q++) {
			num_of_distance_calc[sorted_data[p].getId()]++;
			num_of_distance_calc[sorted_data[q].getId()]++;
			if (reference_point_distance[sorted_data[q].getId()] - reference_point_distance[sorted_data[p].getId()] <= epsilon) {
				e_neighbours[sorted_data[p].getId()].push_back(sorted_data[q].getId());
				e_neighbours[sorted_data[q].getId()].push_back(sorted_data[p].getId());
			} else break;
		}
	}
	// for (int i = 0; i < sample_amount; i++) {
	// 	for (int n : e_neighbours[i]) printf("%i ", n);
	// 	printf("\n");
	// }
	printf("Calculated neighbours for the data using TI\n");	
}

void predecon::calculateVariances() {
	for (int p = 0; p < sample_amount; p++) {
		std::vector<float> vars(attribute_amount, 0);
		for (int q : e_neighbours[p]) {
			for (int a = 0; a < attribute_amount; a++) {
				if (metric == Minkowski2 || metric == Minkowski1){
					vars[a] += pow(data[p].getAttribute(a) - data[q].getAttribute(a), 2)/e_neighbours[p].size();
				}
				else if (metric == MinkowskiInf) {
					vars[a] = std::max((float)pow(data[p].getAttribute(a) - data[q].getAttribute(a), 2)/e_neighbours[p].size(), vars[a]);
				}
			}
		}
		variances.push_back(vars);
	}
	printf("Variances calculated\n");
}

void predecon::calculateSubspacePreferenceVectors() {
	for (int i = 0; i < sample_amount; i++) {
		subspace_preference_vectors.push_back(std::vector<float>(attribute_amount, 0));
		for (int a = 0; a < attribute_amount; a++) {
			subspace_preference_vectors[i][a] = variances[i][a] > delta ? 1 : kappa;
		}
	}
	printf("Calculated subspace preference vectors\n");
}

void predecon::calculatePreferenceWeightedENeighbours() {
	for (int i = 0; i < sample_amount; i++) {
		preference_weighted_e_neighbours.push_back(std::vector<int>(1, i));
	}
	for (int i = 0; i < sample_amount; i++) {
		std::cout << i << '\r';
		calculatePreferenceWeightedENeighbours(i);
	}
	printf("Calculated preference weighted neighbours for the data\n");
}

void predecon::calculatePreferenceWeightedENeighbours(int p) {
	for (int i = p + 1; i < sample_amount; i++) {
		num_of_distance_calc[p]++;
		num_of_distance_calc[p]++;
		num_of_distance_calc[i]++;
		num_of_distance_calc[i]++;
		if (std::max(calculateWeightedDistance(p, i), calculateWeightedDistance(i, p)) <= epsilon) {
			preference_weighted_e_neighbours[p].push_back(i);
			preference_weighted_e_neighbours[i].push_back(p);
		}
	}
}

float predecon::calculateWeightedDistance(int p, int q) {
	float distance = 0;
	if (metric == Minkowski2){
		for (int i = 0; i < attribute_amount; i++) distance += subspace_preference_vectors[p][i] * pow(data[p].getAttribute(i) - data[q].getAttribute(i), 2.0);
		return pow(distance, 0.5);
	} else if (metric == Minkowski1) {
		for (int i = 0; i < attribute_amount; i++) distance += subspace_preference_vectors[p][i] * abs(data[p].getAttribute(i) - data[q].getAttribute(i));
		return distance;
	} else if (metric == MinkowskiInf) {
		for (int i = 0; i < attribute_amount; i++) distance = std::max(subspace_preference_vectors[p][i] * std::abs(data[p].getAttribute(i) - data[q].getAttribute(i)), distance);
		return distance;
	}
}

void predecon::performClustering() {
	current_cluster_id = 0;
	std::vector<int> queue_theta, queue_R;
	for (int i = 0; i < sample_amount; i++) {
		point_type.push_back(0);
		cluster_id.push_back(0);
	}
	// for each sample i
	for (int i = 0; i < sample_amount; i++) {
		// if it's unclassified
		if (cluster_id[i] == 0) {
			// if it's a core
			if (isCore(i)) {
				point_type[i] = 1;
				// generate new id
				current_cluster_id++;
				// add all preference weighted e-neighbours to queue
				for (int n : preference_weighted_e_neighbours[i]) {
					queue_theta.push_back(n);
				}
				// for everything in theta queue
				while (!queue_theta.empty()) {
					// fill R queue
					// 3 DirReach(q,x) conditions (q is first element in theta queue)
					// if q is core
					if (isCore(queue_theta[0])){
						point_type[queue_theta[0]] = 1;
						// and x is in preference weighted e-neighbourhood of q
						for (int n : preference_weighted_e_neighbours[queue_theta[0]]) {
							// and PDIM(x) <= lambda, add x to R queue
							if (calculateSubspacePreferenceDimensionality(n) <= lambda) queue_R.push_back(n);
						}
						// for each element in R
						for (int x : queue_R) {
							// if element is unclassified and not in queue theta put it in there
							if (cluster_id[x] == 0){
								if (std::find(queue_theta.begin(), queue_theta.end(), x) == queue_theta.end()) queue_theta.push_back(x);
								cluster_id[x] = current_cluster_id;
							}
							// if element is unclassified or noise, classify element
							if (point_type[x] == -1) {
								cluster_id[x] = current_cluster_id;
								point_type[x] = 0;
							}
						}
						queue_R.clear();
					}
					queue_theta.erase(queue_theta.begin());
				}
			}
			else {
				point_type[i] = -1;
				cluster_id[i] = -1;
			}
		}
	}
	printf("Clustering complete\n");
}

bool predecon::isCore(int id) {
	return (calculateSubspacePreferenceDimensionality(id) <= lambda && preference_weighted_e_neighbours[id].size() >= mi);
}

int predecon::calculateSubspacePreferenceDimensionality(int i) {
	int dim = 0;
	for (float value : subspace_preference_vectors[i]) if (value == kappa) dim++;
	return dim;
}

void predecon::writeOUTFile() {
	std::ofstream file;
	std::string out_fname("OUT_PREDECON_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_d" + std::to_string(delta) + "_l" + std::to_string(lambda) +
		"_m" + std::to_string(mi) + "_k" + std::to_string(kappa) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + (use_TI ? "_r0" : "") + ".csv");
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

void predecon::writeDEBUGFile() {
	std::ofstream file;
	std::string out_fname("DEBUG_predecon_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_d" + std::to_string(delta) + "_l" + std::to_string(lambda) +
		"_m" + std::to_string(mi) + "_k" + std::to_string(kappa) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + (use_TI ? "_r0" : "") + ".csv");
	file.open(out_fname);
	if (file.is_open()){
		file << "point id, |eps_neighbourhood|, |preference weighted eps_neighbourhood|, PDIM, variance vector, esp_neighbours, preference weighted eps_neighbours\n";
		for (int i = 0; i < sample_amount; i++){
			file << data[i].getId() << ", " << e_neighbours[i].size() << ", " << preference_weighted_e_neighbours[i].size() <<
			", " << calculateSubspacePreferenceDimensionality(i) << ", [";
			for (float v : variances[i]) file << v << (v != variances[i].back() ? "," : "");
			file << "], [";
			for (int n : e_neighbours[i]) file << n << (n != e_neighbours[i].back() ? "," : "");
			file << "], [";
			for (int n : preference_weighted_e_neighbours[i]) file << n << (n != preference_weighted_e_neighbours[i].back() ? "," : "");
			file << "]\n";
		}
		file.close();
	}
	printf("DEBUG file created\n");
}

void predecon::writeSTATFile() {
	std::ofstream file;
	std::string out_fname("STAT_predecon_" + filename + "_D" + std::to_string(attribute_amount) + "_R" +
		std::to_string(sample_amount) + "_e" + std::to_string(epsilon) + "_d" + std::to_string(delta) + "_l" + std::to_string(lambda) +
		"_m" + std::to_string(mi) + "_k" + std::to_string(kappa) + "_mink" +
		((metric == Minkowski1) ? "1" : ((metric == Minkowski2) ? "2" : ((metric == MinkowskiInf) ? "Inf" : "?"))) + (use_TI ? "_r0" : "") + ".csv");
	file.open(out_fname);
	if (file.is_open()){
		file << "Input file: " << filename << "\n";
		file << "Number of dimensions of a point: " << attribute_amount << "\n";
		file << "Number of points: " << sample_amount << "\n";
		file << "Epsilon: " << epsilon << "\n";
		file << "Delta: " << delta << "\n";
		file << "Lambda: " << lambda << "\n";
		file << "Mi: " << mi << "\n";
		file << "Kappa: " << kappa << "\n";
		if (use_TI) {
			file << "Reference point: ";
			for (float v : reference_point) file << v << " ";
				file << "\n";
		}
		file << "Reading a file runtime[us]: " << file_read_time - start_time << "\n";
		if (use_TI) {
			file << "Calculating distance to reference point runtime[us]: " << calculate_reference_point_distance_time - file_read_time << "\n";
			file << "Sorting data due to reference point distance runtime [us]: " << sort_data_time - calculate_reference_point_distance_time << "\n";
			file << "Calculating epsilon neighbours with TI runtime[us]: " << calculate_e_neighbours_time - sort_data_time << "\n";
		} else {
			file << "Calculating epsilon neighbours runtime[us]: " << calculate_e_neighbours_time - file_read_time << "\n";
		}
		file << "Calculating variances runtime[us]: " << calculate_variances_time - calculate_e_neighbours_time << "\n";
		file << "Calculating subspace preference vectors runtime[us]: " << calculate_subspace_preference_vectors_time - calculate_variances_time << "\n";
		file << "Calculating preference weighted epsilon neighbours runtime[us]: " << calculate_preference_weighted_e_neighbours_time - calculate_subspace_preference_vectors_time << "\n";
		file << "Clustering runtime[us]: " << clustering_time - calculate_preference_weighted_e_neighbours_time << "\n";
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
		write_stat_file_time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
		file << "Writing to STAT file runtime[us]: " << write_stat_file_time - write_debug_file_time << "\n";
		file.close();
	}
	printf("STAT file created\n");
}

void predecon::calculateRand() {
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

void predecon::calculatePurity(){
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

void predecon::calculateSilhouetteCoefficient() {
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