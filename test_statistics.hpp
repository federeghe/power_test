#include "distributions.hpp"
#include "simulation.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

inline double get_ad_statistic_upper(const std::vector<double> &F, const std::vector<unsigned int> &f_obs) noexcept {

	double S1=0.0, S2=0.0;

	unsigned int j = 1;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		unsigned int prev_j = j;
		for (; j < prev_j+f_obs[i]; j++) {
			S1 += F[i];
		}
	}
	S1 *= 2.0;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {

		// This is an extreme case, the distribution is completely
		// wrong probably. This may be too pessimistic, but safe.
		if (F[i]==1. && f_obs[i]!=0) return std::numeric_limits<double>::infinity();

		unsigned int prev_j = j;
		for (; j < prev_j+f_obs[i]; j++) {


			double coeff2 = (2. - ((2*j - 1) / config::sample_cardinality) );

			S2 += coeff2 * std::log(1. - F[i]);
		}

	}

	return (double)config::sample_cardinality / 2 - S1 - S2;
}

inline double get_ad_statistic_lower(const std::vector<double> &F, const std::vector<unsigned int> &f_obs) noexcept {

	double S1=0.0, S2=0.0;

	unsigned int j = 1;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {
		unsigned int prev_j = j;
		for (; j < prev_j+f_obs[i]; j++) {
			S1 += F[i];
		}
	}
	S1 *= 2.0;

	for (unsigned int i=0; i < config::sample_cardinality; i++) {

		// This is an extreme case, the distribution is completely
		// wrong probably. This may be too pessimistic, but safe.
		if (F[i]==0. && f_obs[i]!=0) return std::numeric_limits<double>::infinity();

		unsigned int prev_j = j;
		for (; j < prev_j+f_obs[i]; j++) {

			double coeff2 = (2*j - 1) / config::sample_cardinality;

			S2 += coeff2 * std::log(F[i]);
		}
	}

	return - 3 * (double)config::sample_cardinality / 2 + S1 - S2;
}

inline double get_ad_statistic(const std::vector<double> &F, const std::vector<unsigned int> &f_obs) noexcept {
	return get_ad_statistic_upper(F, f_obs) + get_ad_statistic_lower(F, f_obs);
}
