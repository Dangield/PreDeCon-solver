#include <fstream>
#include <cstring>
#include <iostream>

int main(int argc, char const *argv[])
{
	std::ifstream txt_file(argv[1]), pa_file(argv[2]);
	if (!txt_file) {
		printf(".txt file doesn't exist");
		return 1;
	}
	if (!pa_file) {
		printf(".txt file doesn't exist");
		return 1;
	}

	std::ofstream output_file;
	output_file.open("dim512.arff");
	output_file << "@relation dim512\n";
	for (int i = 0; i < 512; i++){
		output_file << "@attribute a" << i << " integer\n";
	}
	output_file << "@attribute class {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}\n";
	output_file << "@data\n";
	
	std::string txt_line, pa_line;
	while (getline(txt_file, txt_line) && getline(pa_file, pa_line)) {
		char prev_c = ' ';
		for (char s : txt_line) {
			if (s == ' ') {
				if (prev_c != ' ') output_file << ",";
			}
			else output_file << s;
			prev_c = s;
		}
		output_file << "," << pa_line << "\n";
	}
}