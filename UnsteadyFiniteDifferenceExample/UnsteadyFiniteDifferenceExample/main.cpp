#include "UnsteadyFiniteDifferenceExample.h"

void test(std::string temporal, std::string index) {
	UnsteadyFiniteDifferenceExample UFDE(1.0, 1.0, 0.1, 40, 0.0001, 1e-8, temporal);
	Data data;
	UFDE.solve(data);
	std::cout << "temporal=" << temporal << " time=" << data.t[data.t.size() - 1] << std::endl;
	data.save(index);
}

void test() {
	test("explicit Euler", "1");
	test("implicit Euler", "2");
	test("Crank-Nicolson", "3");
	test("three-time-level", "4");
}

int main() {
	test();
	return 0;
}
