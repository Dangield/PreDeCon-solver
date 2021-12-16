#include "predecon.hpp"
#include <iostream>
#include <fstream>

predecon::predecon() {
	epsilon = 0;
	delta = 5;
	lambda = 1;
	mi = 1;
	kappa = 100;
	printf("[PREDECON] Predecon solver initialized\n");
}

predecon::~predecon() {
}

void predecon::setParameters(int e, int d, int l, int m, int k = 100) {
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
}

void predecon::printData() {
	printf("\t");
	for (std::string att : attribute_names) {
		printf("%s\t\t", att.c_str());
	}
	printf("\n");

	for (sample single_sample : data) {
		printf("%s\t", single_sample.getId().c_str());
		for (float attribute_value : single_sample.getAttributes()) {
			printf("%3.4f\t\t", attribute_value);
		}
		printf("\n");
	}
}