#include "predecon.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>

predecon::predecon() {
	epsilon = 0;
	delta = 5;
	lambda = 1;
	mi = 1;
	kappa = 100;
	distance_metric = Euclidean;
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
	printf("[PREDECON] Algorythm parameters:\nepsilon = %3.4f\ndelta = %3.4f\nlambda = %d\nmi = %d\nkappa = %3.4f\n", epsilon, delta, lambda, mi, kappa);
}


void predecon::loadDataFromFile(std::string filename){
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
		pos = line.find(delimiter);
		value = line.substr(0, pos);
		data.push_back(sample(value));
		line.erase(0, pos + delimiter.length());
		while ((pos = line.find(delimiter)) != std::string::npos) {
			value = line.substr(0, pos);
			data.back().pushAttribute(std::stof(value));
			line.erase(0, pos + delimiter.length());
		}
		data.back().pushAttribute(std::stof(line));
	}
	attribute_amount = data[0].getAttributeAmount();
}

void predecon::printData() {
	printf("Data\t");
	for (std::string att : attribute_names) {
		printf("%s\t\t", att.c_str());
	}
	if (!e_neighbours.empty()) printf("E-neighbours\t");
	printf("\n");

	for (int i = 0; i < data.size(); i++) {
		printf("%s\t", data[i].getId().c_str());
		for (float attribute_value : data[i].getAttributes()) {
			printf("%3.4f\t\t", attribute_value);
		}
		if (!e_neighbours.empty()) printf("%s\t", std::accumulate(e_neighbours[i].begin(), e_neighbours[i].end(), std::string("")).c_str());
		printf("\n");
	}
}

float predecon::calculateDistance(sample p, sample q) {
	switch(distance_metric){
		case Euclidean:
			return calculateDistanceEuclidean(p, q);
		default:
			return 0;
	}
}

float predecon::calculateDistanceEuclidean(sample p, sample q) {
	float distance = 0;
	for (int i = 0; i < attribute_amount; i++) distance += pow(p.getAttribute(i) - q.getAttribute(i), 2.0);
	// printf("%f", pow(distance, 0.5));
	return pow(distance, 0.5);
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

void predecon::printENeighbours() {
	printf("Data\tE-neighbours\n");
	for (int i = 0; i < data.size(); i++) printf("%s\t%s\n", data[i].getId().c_str(), std::accumulate(e_neighbours[i].begin(), e_neighbours[i].end(), std::string("")).c_str());
}