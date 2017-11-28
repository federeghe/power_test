#include <cmath>
#include <cassert>


#include "montecarlo.hpp"
#include "simulation.hpp"

static inline void put_output(double x, std::vector<unsigned int> &output) {

		x = std::ceil(x * (1.0/config::step)) / (1.0/config::step);
		// Get the corresponding index
		double index = (x - config::dist_min) / config::step;

		if (index < 0) index = 0;
		if (index > config::size-1) index = config::size-1;

		output[index]++;
}

extern void calculate_cumulative(const std::vector<unsigned int> &input, std::vector<double> &output) {
	unsigned int total = 0;
	for (unsigned int i=0; i<input.size(); i++) {
		total += input[i];
		output[i] = ((double)total) / config::sample_cardinality; 
	}
}


std::mt19937 rng;

template <>
void montecarlo_frequencies<NORMAL>(std::vector<unsigned int> &output) {

	assert(output.size() >= config::sample_cardinality);
	std::normal_distribution<double> distribution(0.0,1.0);

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = distribution(rng);
		put_output(x, output);
	}

}

template <>
void montecarlo_frequencies<UNIFORM>(std::vector<unsigned int> &output) {

	assert(output.size() >= config::sample_cardinality);
	std::uniform_real_distribution<double> distribution(config::dist_min, config::dist_max-1);

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = distribution(rng);
		put_output(x, output);
	}

}


