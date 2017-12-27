#ifndef MONTECARLO_HPP_
#define MONTECARLO_HPP_

#include "simulation.hpp"

#include <vector>
#include <random>

typedef enum distribution_e {
	NORMAL,
	UNIFORM,
	T_STUDENT,
	GEV

} distribution_t;


inline void put_output(double x, std::vector<unsigned int> &output) {
	double step = (config::dist_max - config::dist_min) / config::sample_cardinality;
	double ratio = 1.0/step;

	x = std::ceil(x * ratio) / ratio;
	// Get the corresponding index
	double index = (x - config::dist_min) / step;

	if (index < 0) index = 0;
	if (index > config::sample_cardinality-1) index = config::sample_cardinality-1;

	output[index]++;
}

template <distribution_t T=NORMAL>
extern void montecarlo_frequencies(std::mt19937 &rng, std::vector<unsigned int> &output);

template <int num=0, int den=1, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline double montecarlo_evt_sample(std::mt19937 &rng) {

	static_assert(den != 0, "Denominator must be != from 0!");

	constexpr double xi=((double)num)/((double)den);
	constexpr double mu=((double)mu_num)/((double)mu_den);
	constexpr double sigma=((double)sigma_num)/((double)sigma_den);

	std::uniform_real_distribution<double> distribution(mu,sigma);

	double x = distribution(rng);

	// Convert it with the p-quantile
	if (num != 0) {
		x = mu + sigma * (1-std::pow((-std::log(x)),(-xi)))/(-xi);
	} else {
		x = mu - sigma * std::log(-std::log(x));
	}

	return x;

}


template <int num=0, int den=1, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline void montecarlo_frequencies_evt(std::mt19937 &rng, std::vector<unsigned int> &output) {
	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double x = montecarlo_evt_sample<num, den, mu_num, mu_den, sigma_num, sigma_den>(rng);
		put_output(x, output);
	}
}

extern void calculate_cumulative(const std::vector<unsigned int> &input, std::vector<double> &output);

extern std::mt19937 rng;

#endif
