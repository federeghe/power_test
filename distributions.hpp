#ifndef DISTRIBUTIONS_HPP_
#define DISTRIBUTIONS_HPP_

#include <cmath>
#include <vector>

extern void fill_evt_cumulative(double xi, double mu, double sigma, std::vector<double> &output);
extern void fill_evt_cumulative_xi0(double mu, double sigma, std::vector<double> &output);

template <int num, int den, int mu_num, int mu_den, int sigma_num, int sigma_den>
extern double get_evt_cumulative(double x) {

	static_assert(den != 0, "Denominator must be != from 0!");
	static_assert(mu_den != 0, "Denominator must be != from 0!");
	static_assert(sigma_den != 0, "Denominator must be != from 0!");

	constexpr double xi=((double)num)/((double)den);
	constexpr double mu=((double)mu_num)/((double)mu_den);
	constexpr double sigma=((double)sigma_num)/((double)sigma_den);

	double t_x;

	if (num != 0) {
		double cond_value =  1.0 + xi * ((x-mu)/sigma);

		if ( cond_value < 0 ) {
			return xi > 0 ? 0 : 1;
		}

		// t(x) calculation
		t_x = std::pow(cond_value,  - 1.0 / xi);
	} else {
		// t(x) calculation
		t_x = std::exp( - ( x - mu ) / sigma);
	}

	// CDF calculation
	return std::exp(-t_x);
}

#endif
