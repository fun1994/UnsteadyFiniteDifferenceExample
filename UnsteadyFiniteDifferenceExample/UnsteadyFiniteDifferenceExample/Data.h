#pragma once
#include <iostream>
#include <vector>
#include <fstream>

class Data {
	void save(std::vector<double>& data, std::string filename) {
		std::ofstream file("./data/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			file << data[i];
			if (i < data.size() - 1) {
				file << " ";
			}
		}
		file.close();
	}
public:
	std::vector<double> x;
	std::vector<double> t;
	std::vector<std::vector<double>> phi;
	void save(std::string index) {
		save(x, "x");
		save(phi[phi.size() - 1], "phi" + index);
	}
};
