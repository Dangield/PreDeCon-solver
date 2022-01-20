#include <vector>
#include <string>

#ifndef SAMPLE_H
#define SAMPLE_H

class sample
{
private:
	int id;
	std::vector<float> attribute_values;
	static int global_id;
public:
	sample();
	~sample();
	void pushAttribute(float);
	int getId();
	std::vector<float> getAttributes();
	float getAttribute(int);
	int getAttributeAmount();
};

#endif