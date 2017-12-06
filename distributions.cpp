#include <cmath>
#include <cassert>

#include "distributions.hpp"
#include "simulation.hpp"

extern void fill_evt_cumulative(double xi, double mu, double sigma, std::vector<double> &output) {

	assert(output.size() >= config::size);

	for (unsigned int i=0; i < config::size; i++) {

		double x = config::dist_min + i * config::step;
		double cond_value =  1.0 + xi * ((x-mu)/sigma);

		if ( cond_value < 0 ) {
			output[i] = 0;
			continue;
		}

		// t(x) calculation
		double t_x = std::pow(cond_value,  - 1.0 / xi);

		// CDF calculation
		output[i] = std::exp(-t_x);
	}
}

extern void fill_evt_cumulative_xi0(double mu, double sigma, std::vector<double> &output) {

	assert(output.size() >= config::size);

	for (unsigned int i=0; i < config::size; i++) {

		double x = config::dist_min + i * config::step;

		// t(x) calculation
		double t_x = std::exp( - ( x - mu ) / sigma);

		// CDF calculation
		output[i] = std::exp(-t_x);
	}

}
