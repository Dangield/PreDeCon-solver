#include "sample.hpp"

int sample::global_id = 0;

sample::sample(){
	id = global_id;
	global_id++;
}

sample::sample(std::string id_n) {
	id = id_n;
	global_id++;
}

sample::~sample(){
}

void sample::pushAttribute(float value) {
	attribute_values.push_back(value);
}

std::string sample::getId() {
	return id;
}

std::vector<float> sample::getAttributes() {
	return attribute_values;
}

float sample::getAttribute(int n){
	return attribute_values[n];
}

int sample::getAttributeAmount() {
	return attribute_values.size();
}