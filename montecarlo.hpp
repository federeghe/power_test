#ifndef MONTECARLO_HPP_
#define MONTECARLO_HPP_

#include <vector>
#include <random>

typedef enum distribution_e {
	NORMAL,
	UNIFORM,
	EVT_I,
	EVT_II,
	EVT_III,

} distribution_t;

template <distribution_t T=NORMAL>
extern void montecarlo_frequencies(std::vector<unsigned int> &output);

extern void calculate_cumulative(const std::vector<unsigned int> &input, std::vector<double> &output);

extern std::mt19937 rng;

#endif
