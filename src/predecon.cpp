#include "predecon.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cstring>
#include <algorithm>
#include "matplotlibcpp.h"

predecon::predecon() {
	epsilon = 0;
	delta = 5;
	lambda = 1;
	mi = 1;
	kappa = 100;
	distance_metric = Euclidean;
	current_cluster_id = 0;
	printf("[PREDECON] Predecon solver initialized\n");
}

predecon::~predecon() {
}

void predecon::setParameters(float e, float d, int l, int m, float k) {
	epsilon = e;
	delta = d;
	lambda = l;
	mi = m;
	kappa = k;
}

void predecon::printParameters() {
	printf("[PREDECON] Algorythm parameters:\nepsilon = %.3f\ndelta = %.3f\nlambda = %d\nmi = %d\nkappa = %.3f\n", epsilon, delta, lambda, mi, kappa);
}


void predecon::loadDataFromFile(std::string filename, bool dataHasIds){
	// open data file
	std::ifstream file(filename);
	if (!file) {
		printf("[PREDECON] Data file doesn't exist!\n");
		return;
	}
	std::string line;
	std::string delimiter = ";";
	std::string value;
	size_t pos = 0;
	// read attribute names
	getline(file, line);
	while ((pos = line.find(delimiter)) != std::string::npos) {
		value = line.substr(0, pos);
		attribute_names.push_back(value);
		line.erase(0, pos + delimiter.length());
	}
	attribute_names.push_back(line);
	// for (std::string att : attribute_names) std::cout << att << std::endl;
	// read samples
	while (getline(file, line)) {
		// std::cout << line << std::endl;
		if (dataHasIds){
			pos = line.find(delimiter);
			value = line.substr(0, pos);
			data.push_back(sample(value));
			line.erase(0, pos + delimiter.length());
		}
		else data.push_back(sample());
		while ((pos = line.find(delimiter)) != std::string::npos) {
			value = line.substr(0, pos);
			data.back().pushAttribute(std::stof(value));
			line.erase(0, pos + delimiter.length());
		}
		data.back().pushAttribute(std::stof(line));
	}
	attribute_amount = data[0].getAttributeAmount();
}

void predecon::printData(bool waitForInput) {
	std::cout << "\033[2J\033[1;1H";
	printf("Data\t");
	for (std::string att : attribute_names) {
		printf("%s\t\t", att.c_str());
	}
	if (!e_neighbours.empty()) printf("E-neighbours\t");
	if (!variances.empty()) printf("Variances\t");
	if (!subspace_preference_vectors.empty()) printf("SPVectors\t");
	if (!preference_weighted_e_neighbours.empty()) printf("PWE-neighbours\t");
	if (!core.empty()) printf("Core\t");
	if (!cluster.empty()) printf("ClusterId\t");
	printf("\n");

	for (int i = 0; i < data.size(); i++) {
		printf("%s\t", data[i].getId().c_str());
		for (float attribute_value : data[i].getAttributes()) {
			printf("%.3f\t\t", attribute_value);
		}
		if (!e_neighbours.empty()) printf("%s\t\t", std::accumulate(e_neighbours[i].begin(), e_neighbours[i].end(), std::string("")).c_str());
		if (!variances.empty()) for (float var : variances[i]) printf("%.3f\t", var);
		if (!subspace_preference_vectors.empty()) for (float w : subspace_preference_vectors[i]) printf("%.3f\t", w);
		if (!preference_weighted_e_neighbours.empty()) printf("%s\t\t", std::accumulate(preference_weighted_e_neighbours[i].begin(), preference_weighted_e_neighbours[i].end(), std::string("")).c_str());
		if (!core.empty()) printf("%s\t", property_strings[core[i]].c_str());
		if (!cluster.empty()) if (cluster[i] > 0) printf("%d\t", cluster[i]); else printf("-\t");
		printf("\n");
	}

	if (!queue_theta.empty()) {
		printf("Queue theta:\t");
		for (int id : queue_theta) printf("%s ", data[id].getId().c_str());
		printf("\n");
	}

	if (!queue_R.empty()) {
		printf("Queue R:\t");
		for (int id : queue_R) printf("%s ", data[id].getId().c_str());
		printf("\n");
	}

	if (waitForInput) std::cin.get();
}

float predecon::calculateDistance(sample p, sample q) {
	switch(distance_metric){
		case Euclidean:
			return calculateDistanceEuclidean(p, q);
		case Minkowsky1:
			return calculateDistanceMinkowsky1(p, q);
		case Minkowsky2:
			return calculateDistanceMinkowsky2(p, q);
		case MinkowskyInf:
			return calculateDistanceMinkowsky2(p, q);
		default:
			return 0;
	}
}

float predecon::calculateDistanceEuclidean(sample p, sample q) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += pow(p.getAttribute(i) - q.getAttribute(i), 2.0);
	return pow(distance, 0.5);
}

float predecon::calculateDistanceMinkowsky1(sample p, sample q) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += abs(p.getAttribute(i) - q.getAttribute(i));
	return distance;
}

float predecon::calculateDistanceMinkowsky2(sample p, sample q) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += pow(abs(p.getAttribute(i) - q.getAttribute(i)), 2.0);
	return pow(distance, 0.5);
}

float predecon::calculateDistanceMinkowskyInf(sample p, sample q) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance = std::max(float(abs(p.getAttribute(i) - q.getAttribute(i))), distance);
	return distance;
}

void predecon::setDistanceMetric(DistanceMetric m) {
	distance_metric = m;
}

void predecon::calculateENeighbours() {
	for (sample p : data) e_neighbours.push_back(calculateENeighbours(p));
}

std::vector<std::string> predecon::calculateENeighbours(sample p) {
	std::vector<std::string> neighbours;
	for (sample q : data) if (calculateDistance(p, q) <= epsilon) neighbours.push_back(q.getId());
	return neighbours;
}

void predecon::calculateVariances() {
	for (int i = 0; i < data.size(); i++) {
		variances.push_back(calculateVariances(data[i], e_neighbours[i]));
	}
	// e_neighbours.clear();
}

std::vector<float> predecon::calculateVariances(sample p, std::vector<std::string> neighbours) {
	std::vector<float> vars;
	for (int i = 0; i < attribute_amount; i++) vars.push_back(0.0);
	for (std::string neighbour : neighbours) {
		for (sample q : data) {
			if (strcmp(q.getId().c_str(), neighbour.c_str()) == 0) {
				for (int i = 0; i < attribute_amount; i++) {
					vars[i] += pow(p.getAttribute(i) - q.getAttribute(i), 2.0);
				}
				break;
			}
		}
	}
	for (int i = 0; i < attribute_amount; i++) vars[i] = vars[i]/neighbours.size();
	return vars;
}

void predecon::calculateSubspacePreferenceVectors() {
	for (std::vector<float> vars : variances) subspace_preference_vectors.push_back(calculateSubspacePreferenceVectors(vars));
	// variances.clear();
}

std::vector<float> predecon::calculateSubspacePreferenceVectors(std::vector<float>  vars) {
	std::vector<float> v;
	for (float var : vars) if (var > delta) v.push_back(1); else v.push_back(kappa);
	return v;
}

void predecon::calculatePreferenceWeightedENeighbours() {
	for (int i = 0; i < data.size(); i++) preference_weighted_e_neighbours.push_back(calculatePreferenceWeightedENeighbours(data[i], subspace_preference_vectors[i]));
}

std::vector<std::string> predecon::calculatePreferenceWeightedENeighbours(sample p, std::vector<float> w) {
	std::vector<std::string> neighbours;
	for (int i = 0; i < data.size(); i++) if (std::max(calculateWeightedDistance(p, data[i], w), calculateWeightedDistance(data[i], p, subspace_preference_vectors[i])) <= epsilon) neighbours.push_back(data[i].getId());
	return neighbours;
}

float predecon::calculateWeightedDistance(sample p, sample q, std::vector<float> w) {
	switch(distance_metric){
		case Euclidean:
			return calculateWeightedDistanceEuclidean(p, q, w);
		case Minkowsky1:
			return calculateWeightedDistanceMinkowsky1(p, q, w);
		case Minkowsky2:
			return calculateWeightedDistanceMinkowsky2(p, q, w);
		case MinkowskyInf:
			return calculateWeightedDistanceMinkowskyInf(p, q, w);
		default:
			return 0;
	}
}

float predecon::calculateWeightedDistanceEuclidean(sample p, sample q, std::vector<float> w) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += w[i] * pow(p.getAttribute(i) - q.getAttribute(i), 2.0);
	return pow(distance, 0.5);
}

float predecon::calculateWeightedDistanceMinkowsky1(sample p, sample q, std::vector<float> w) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += w[i] * abs(p.getAttribute(i) - q.getAttribute(i));
	return distance;
}

float predecon::calculateWeightedDistanceMinkowsky2(sample p, sample q, std::vector<float> w) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += w[i] * pow(abs(p.getAttribute(i) - q.getAttribute(i)), 2.0);
	return pow(distance, 0.5);
}

float predecon::calculateWeightedDistanceMinkowskyInf(sample p, sample q, std::vector<float> w) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance = std::max(w[i] * abs(p.getAttribute(i) - q.getAttribute(i)), distance);
	return distance;
}

void predecon::solve() {
	printData(true);
	for (int i = 0; i < data.size(); i++) {
		core.push_back(Unknown);
		cluster.push_back(0);
	}
	// for each sample i
	for (int i = 0; i < data.size(); i++) {
		printData();
		printf("Sample %s\n", data[i].getId().c_str());
		std::cin.get();
		// if it's unclassified
		if (cluster[i] == 0) {
			// if it's a core
			if (isCore(i)) {
				// generate new id
				current_cluster_id++;
				// add all preference weighted e-neighbours to queue
				for (std::string neighbour : preference_weighted_e_neighbours[i]) {
					for (int j = 0; j < data.size(); j++) {
						if (strcmp(data[j].getId().c_str(), neighbour.c_str()) == 0) {
							queue_theta.push_back(j);
						}
					}
				}
				printData(true);
				// for everything in theta queue
				while (!queue_theta.empty()) {
					// fill R queue
					// 3 DirReach(q,x) conditions (q is first element in theta queue)
					// if q is core
					if (isCore(queue_theta[0])){
						// and x is in preference weighted e-neighbourhood of q
						for (std::string neighbour : preference_weighted_e_neighbours[queue_theta[0]]) {
							for (int j = 0; j < data.size(); j++){
								if (strcmp(data[j].getId().c_str(), neighbour.c_str()) == 0) {
									// and PDIM(x) <= lambda, add x to R queue
									if (calculateSubspacePreferenceDimensionality(j) <= lambda) queue_R.push_back(j);
								}
							}
						}
						printData(true);
						// for each element in R
						for (int id : queue_R) {
							// if element is unclassified and not in queue theta put it in there
							if (cluster[id] == 0 && std::find(queue_theta.begin(), queue_theta.end(), id) == queue_theta.end()) queue_theta.push_back(id);
							// if element is unclassified or noise, classify element
							if (cluster[id] == 0 || isNoise(id)) cluster[id] = current_cluster_id;
						}
						queue_R.clear();
					}
					queue_theta.erase(queue_theta.begin());
					printData(true);
				}
			}
			else core[i] = Noise;
		}
	}
}

bool predecon::isCore(int id) {
	switch(core[id]){
		case Unknown:
			if (calculateSubspacePreferenceDimensionality(id) <= lambda && preference_weighted_e_neighbours[id].size() >= mi){
				core[id] = Core;
				return true;
			}
		case Noise:
			return false;
		case Core:
			return true;
	}
}

bool predecon::isNoise(int id) {
	switch(core[id]){
		case Unknown:
			if (calculateSubspacePreferenceDimensionality(id) <= lambda && preference_weighted_e_neighbours[id].size() >= mi) {
				core[id] = Core;
				return false;
			}
		case Noise:
			return true;
		case Core:
			return false;
	}
}

int predecon::calculateSubspacePreferenceDimensionality(int id) {
	int dim = 0;
	for (float value : subspace_preference_vectors[id]) if (value == kappa) dim++;
	return dim;
}

void predecon::plotData() {
	if (attribute_amount == 2){
		std::vector<float> x, y;
		for (int i = 0; i < data.size(); i++) {
			x.push_back(data[i].getAttribute(0));
			y.push_back(data[i].getAttribute(1));
		}
		matplotlibcpp::plot(x, y, "*");
		if (!subspace_preference_vectors.empty() && distance_metric == Euclidean || distance_metric == Minkowsky2) for (int i = 0; i < data.size(); i++) {
			std::vector<float> elipse_x, elipse_y;
			for (float angle = 0; angle <= 2*M_PI; angle += 2*M_PI/360) {
				elipse_x.push_back(cos(angle)*epsilon/pow(subspace_preference_vectors[i][0], 0.5) + x[i]);
				elipse_y.push_back(sin(angle)*epsilon/pow(subspace_preference_vectors[i][1], 0.5) + y[i]);
			}
			matplotlibcpp::plot(elipse_x, elipse_y, "r");
		}
		if (!cluster.empty()) for (int i = 0; i < cluster.size(); i++) matplotlibcpp::text(x[i], y[i], data[i].getId() + "|" + std::to_string(cluster[i]));
		matplotlibcpp::axis("equal");
		matplotlibcpp::show();
	}
}