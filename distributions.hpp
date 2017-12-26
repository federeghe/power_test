#ifndef DISTRIBUTIONS_HPP_
#define DISTRIBUTIONS_HPP_

#include <cmath>
#include <vector>

extern void fill_evt_cumulative(double xi, double mu, double sigma, std::vector<double> &output);
extern void fill_evt_cumulative_xi0(double mu, double sigma, std::vector<double> &output);

template <int XI_NUM, int XI_DEN>
inline double get_evt(double mu, double sigma, double x) noexcept {

	constexpr double xi = ((double)XI_NUM) / ((double)XI_DEN);
	double t_x;

	if constexpr (XI_NUM != 0.) {
		// t(x) calculation
		double cond_value =  1.0 + xi * ((x-mu)/sigma);
		if (cond_value < 0) {
			return xi > 0 ? 0 : 1;
		}
		t_x = std::pow(cond_value,  - 1.0 / xi);
		
	} else {
		// t(x) calculation
		t_x = std::exp( - ( x - mu ) / sigma);
	}

	// CDF calculation
	return std::exp(-t_x);
}

#endif
