#include <vector>
#include <string>

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