#include "distributions.hpp"
#include "simulation.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

template <int num, int den, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline double get_ad_statistic_upper(const double *samples) noexcept {

	double S1=0.0, S2=0.0;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double coeff1 = (2. - ((2.*((int)i+1.) - 1.) / ((double)config::sample_cardinality)) );
		double Fi = get_evt_cumulative<num, den, mu_num, mu_den, sigma_num, sigma_den>(samples[i]);

		// This is an extreme case, the distribution is completely
		// wrong probably. This may be too pessimistic, but safe.
		if (Fi == 1.) return std::numeric_limits<double>::infinity();

		S1 += Fi;
		S2 += coeff1 * std::log(1. - Fi);

	}

	S1 *= 2.0;

	return (double)config::sample_cardinality / 2 - S1 - S2;
}

template <int num, int den, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline double get_ad_statistic_lower(const double *samples) noexcept {

	double S1=0.0, S2=0.0;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		double coeff2 = (((2.*((int)i+1) - 1.) / ((double)config::sample_cardinality)) );
		double Fi = get_evt_cumulative<num, den, mu_num, mu_den, sigma_num, sigma_den>(samples[i]);

		// This is an extreme case, the distribution is completely
		// wrong probably. This may be too pessimistic, but safe.
		if (Fi == 0.) return std::numeric_limits<double>::infinity();

		S1 += Fi;
		S2 += coeff2 * std::log(Fi);

	}

	S1 *= 2.0;

	return - 3 * (double)config::sample_cardinality / 2 + S1 - S2;
}

template <int num, int den, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline double get_ad_statistic(const double *samples) noexcept {
	double upper = get_ad_statistic_upper<num, den, mu_num, mu_den, sigma_num, sigma_den>(samples);
	double lower = get_ad_statistic_lower<num, den, mu_num, mu_den, sigma_num, sigma_den>(samples);

	return upper + lower;
}



template <int num, int den, int mu_num=0, int mu_den=1, int sigma_num=1, int sigma_den=1>
inline double get_ad_statistic_full(const double *samples) noexcept {

	double n = config::sample_cardinality;

	double S1 = 0.0;
	double S2 = 0.0;

	for (unsigned int i=0; i < n; i++) {

		int coeff1 = 2*(((int)i)+1)-1;
		int coeff2 = 2*n - 2*(((int)i)+1) + 1;
		double Fi = get_evt_cumulative<num, den, mu_num, mu_den, sigma_num, sigma_den>(samples[i]);

		S1 += coeff1 * std::log(Fi);
		S2 += coeff2 * std::log(1.-Fi);
	}


	return - n - (1./n * S1) - (1./n * S2);
}




