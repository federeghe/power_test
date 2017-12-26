#include <cmath>
#include <cassert>


#include "montecarlo.hpp"


extern void calculate_cumulative(const std::vector<unsigned int> &input, std::vector<double> &output) {
	unsigned int total = 0;
	for (unsigned int i=0; i<input.size(); i++) {
		total += input[i];
		output[i] = ((double)total) / config::sample_cardinality; 
	}
}

template <>
void montecarlo_frequencies<NORMAL>(std::mt19937 &rng, std::vector<unsigned int> &output) {

	assert(output.size() >= config::sample_cardinality);
	std::normal_distribution<double> distribution(0.0,1.0);

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = distribution(rng);
		put_output(x, output);
	}

}

template <>
void montecarlo_frequencies<UNIFORM>(std::mt19937 &rng, std::vector<unsigned int> &output) {

	assert(output.size() >= config::sample_cardinality);
	std::uniform_real_distribution<double> distribution(-2,3);

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = distribution(rng);
		put_output(x, output);
	}

}

template <>
void montecarlo_frequencies<T_STUDENT>(std::mt19937 &rng, std::vector<unsigned int> &output) {

	assert(output.size() >= config::sample_cardinality);
	std::student_t_distribution<double> distribution(10.0);

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = distribution(rng);
		put_output(x, output);
	}

}

