#include "test_statistics.hpp"
#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <omp.h>

#define NR_RUNS 1000
#define NR_RUNS_CRIT 1000
#define DEBUG 0

#define TEST_KS  1
#define TEST_AD  2 
#define TEST_MAD 3

#define MAX_SAMPLE_CARDINALITY 20000

unsigned int config::sample_cardinality;

thread_local double sample[MAX_SAMPLE_CARDINALITY];

template <int num, int den, int significance_num, int significance_den>
double compute_ad_crit_value() {
	double crits[NR_RUNS_CRIT];

	#pragma omp parallel
	{
		thread_local static std::mt19937 random_gen(std::random_device{}());
		#pragma omp for
		for (unsigned long i=0; i < NR_RUNS_CRIT; i++) {
			for (unsigned int j=0; j < config::sample_cardinality; j++) {
				sample[j] = montecarlo_evt_sample<num, den>(random_gen);
			}
			std::sort(std::begin(sample), std::end(sample));
#if TEST==TEST_AD
			crits[i] = get_ad_statistic<num, den>(sample);
#endif
#if TEST==TEST_MAD
			crits[i] = get_ad_statistic_upper<num, den>(sample);
#endif
		}
	}
	std::sort(std::begin(crits),std::end(crits));

	return crits[(int)(((double)significance_num/(double)significance_den) * NR_RUNS_CRIT)];
}

#ifndef TEST
	#error TEST must be defined
#endif

#if TEST==TEST_AD || TEST==TEST_MAD
template <int num, int den>
static inline void run() noexcept {
	constexpr double xi = ((double)num) / ((double)den);

	const double ad_critical_value1 = compute_ad_crit_value<0, 1, 99, 100>();
	const double ad_critical_value2 = compute_ad_crit_value<0, 1, 95, 100>();

	std::cout << config::sample_cardinality << " " << xi << " ";

	unsigned long reject1=0;
	unsigned long reject2=0;

	#pragma omp parallel firstprivate(ad_critical_value1, ad_critical_value2, config::sample_cardinality)
	{
		thread_local static std::mt19937 random_gen(std::random_device{}());

		#pragma omp for reduction(+:reject1,reject2)
		for (unsigned long j=0; j < NR_RUNS; j++) {


			for (unsigned int i=0; i < config::sample_cardinality; i++) {
				sample[i] = montecarlo_evt_sample<num, den>(random_gen);
			}

			std::sort(std::begin(sample), std::end(sample));

#if TEST == TEST_AD
			double S = get_ad_statistic<0, 1>(sample);
#elif TEST == TEST_MAD
			double S = get_ad_statistic_upper<0, 1>(sample);
#endif

			if ( S > ad_critical_value1 ) {
				reject1++;
			}
			if ( S > ad_critical_value2 ) {
				reject2++;
			}

		}
	
	}
	std::cout << reject1;
	std::cout << " ";
	std::cout << reject2;
	std::cout << std::endl;

}
#endif

#if TEST==1
template <int num, int den>
static inline void run() noexcept {
	constexpr double xi = ((double)num) / ((double)den);

	std::vector<double> F_evt;

	F_evt.resize(config::sample_cardinality);

	fill_evt_cumulative_xi0(0, 1, F_evt);


	double ks_critical_value1 = get_ks_critical<1, 100>(config::sample_cardinality);
	double ks_critical_value2 = get_ks_critical<5, 100>(config::sample_cardinality);

	std::cout << config::sample_cardinality << " " << xi << " ";

	unsigned long reject1=0;
	unsigned long reject2=0;

	#pragma omp parallel firstprivate(F_evt)
	{
		thread_local static std::mt19937 random_gen(std::random_device{}());


		std::vector<unsigned int> freq_montecarlo;
		std::vector<double> F_montecarlo;
		freq_montecarlo.resize(config::sample_cardinality);
		F_montecarlo.resize(config::sample_cardinality);

		#pragma omp for reduction(+:reject1,reject2)
		for (unsigned long j=0; j < NR_RUNS; j++) {

			std::fill(freq_montecarlo.begin(), freq_montecarlo.end(), 0);

			montecarlo_frequencies_evt<num,den>(random_gen, freq_montecarlo);

			calculate_cumulative(freq_montecarlo, F_montecarlo);


			bool b_reject1=false;
			bool b_reject2=false;
			for (unsigned int i=0; i < config::sample_cardinality; i++) {
				double diff = std::abs(F_evt[i] - F_montecarlo[i]);
//				std::cout << diff << ">" << ks_critical_value1 << std::endl;
				if (diff > ks_critical_value1) {
					b_reject1 = true;
//					std::cout << "    RIGETTONE 1 - " << F_evt[i-1] <<"==" << F_montecarlo[i-1] << " (" << freq_montecarlo[i-1] << ")" << std::endl;
//					std::cout << "    RIGETTONE 1 - " << F_evt[i] <<"==" << F_montecarlo[i] << " (" << freq_montecarlo[i] << ")"<< std::endl;
//					std::cout << "    RIGETTONE 1 - " << F_evt[i+1] <<"==" << F_montecarlo[i+1] << " (" << freq_montecarlo[i+1] << ")"<< std::endl;
				}
				if (diff > ks_critical_value2) {
					b_reject2 = true;
//					std::cout << "    RIGETTONE 2 - i=" << i << " " << F_evt[i-2] <<"==" << F_montecarlo[i-1] << " (" << freq_montecarlo[i-2] << ")" << std::endl;
//					std::cout << "    RIGETTONE 2 - i=" << i << " " << F_evt[i-1] <<"==" << F_montecarlo[i-1] << " (" << freq_montecarlo[i-1] << ")" << std::endl;
//					std::cout << "    RIGETTONE 2 - i=" << i << " " << F_evt[i] <<"==" << F_montecarlo[i] << " (" << freq_montecarlo[i] << ")"<< std::endl;
//					std::cout << "    RIGETTONE 2 - i=" << i << " " << F_evt[i+1] <<"==" << F_montecarlo[i+1] << " (" << freq_montecarlo[i+1] << ")"<< std::endl;
//					std::cout << "    RIGETTONE 2 - i=" << i << " " << F_evt[i+2] <<"==" << F_montecarlo[i+2] << " (" << freq_montecarlo[i+2] << ")" << std::endl;
				}
				if (diff > ks_critical_value1 && diff > ks_critical_value2) break;
			}
			if (b_reject1) reject1++;
			if (b_reject2) reject2++;


		}
	
	}
	std::cout << reject1;
	std::cout << " ";
	std::cout << reject2;
	std::cout << std::endl;

}
#endif

template<int x, int to>
struct static_for
{
	void operator()() {
		run<x, 1000>();
		static_for<x+1,to>()();
	}
};

template<int to>
struct static_for<to,to>
{
    void operator()() 
    {}
};


int main(int argc, char* argv[]) {



	std::vector<unsigned int> cardinalities;
	for (int i=0; i<500; i++) {
		cardinalities.push_back((i+1)*10);
	}


	for (const auto i : cardinalities) {
		config::sample_cardinality = i;
		assert(config::sample_cardinality <= MAX_SAMPLE_CARDINALITY);
		static_for<1,501>()();
	}
	
	return 0;
}




