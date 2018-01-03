#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"
#include "test_statistics.hpp"

#include <algorithm>
#include <iomanip>    
#include <iostream>
#include <omp.h>

#define DIST_T_STUDENT 1
#define DIST_UNIFORM 2
#define DIST_EVT0 3
#define DIST_EVT_0_5 4
#define DIST_EVT_N_0_5 5
#define DIST_NORMAL 6

//#define H1 DIST_EVT_0_5

//#define H0 DIST_EVT_0_5

//#define DEBUG 1

#define TEST_KS 1
#define TEST_AD 2

#if H0 == DIST_EVT0
	#define EVT_PARAM_NUM 0
	#define EVT_PARAM_DEN 1
constexpr double evt_param = 0.0;
#elif H0 == DIST_EVT_0_5
	#define EVT_PARAM_NUM 1
	#define EVT_PARAM_DEN 2
constexpr double evt_param = 0.5;
#elif H0 == DIST_EVT_N_0_5
	#define EVT_PARAM_NUM -1
	#define EVT_PARAM_DEN 2
constexpr double evt_param = -0.5;
#else
	#error "Unknown distribution"
#endif

#if H1 == DIST_T_STUDENT
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<T_STUDENT>
	#define MONTE_CARLO_SAMPLING montecarlo_t_student_sample
#elif H1 == DIST_UNIFORM
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<UNIFORM>
	#define MONTE_CARLO_SAMPLING montecarlo_uniform_sample
#elif H1 == DIST_NORMAL
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<NORMAL>
	#define MONTE_CARLO_SAMPLING montecarlo_normal_sample
#elif H1 == DIST_EVT0
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<0,1>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<0, 1>
#elif H1 == DIST_EVT_0_5
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<1,2>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<1, 2>
#elif H1 == DIST_EVT_N_0_5
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<-1,2>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<-1, 2>
#else
	#error "Unknown"
#endif

unsigned int config::sample_cardinality;

template<int alpha_num, int alpha_den>
unsigned long perform_ad_run() {

	double ad_critical_value = AndersonDarlingCV_EV<alpha_num, alpha_den>::get_critical_value(config::sample_cardinality, evt_param);

	unsigned long reject=0;

	#pragma omp parallel
	{
		std::mt19937 random_gen(std::random_device{}());
		std::vector<double> sample;
		sample.resize(config::sample_cardinality);

		#pragma omp for reduction(+:reject)
		for (unsigned long j=0; j < config::runs; j++) {

			for (unsigned int i=0; i < config::sample_cardinality; i++) {
				sample[i] = MONTE_CARLO_SAMPLING(random_gen);
			}

			std::sort(sample.begin(), sample.end());

			double S = get_ad_statistic<EVT_PARAM_NUM, EVT_PARAM_DEN>(sample);
			if ( S > ad_critical_value ) {
				reject++;
			}
#if DEBUG
			std::cout << "S=" << S << " crit_value=" << ad_critical_value << std::endl;
#endif
		}
	}
	return reject;
}

template<int alpha_num, int alpha_den>
unsigned long perform_ks_run() {
	std::vector<double> F_evt;

	F_evt.resize(config::sample_cardinality);

	#if H0 == DIST_EVT0
	fill_evt_cumulative_xi0(0, 1, F_evt);
	#else
	fill_evt_cumulative(evt_param, 0, 1, F_evt);
	#endif

	double ks_critical_value = get_ks_critical<alpha_num, alpha_den>(config::sample_cardinality);

	unsigned long reject=0;

	#pragma omp parallel firstprivate(F_evt)
	{

		thread_local static std::mt19937 random_gen(std::random_device{}());


		std::vector<unsigned int> freq_montecarlo;
		std::vector<double> F_montecarlo;
		freq_montecarlo.resize(config::sample_cardinality);
		F_montecarlo.resize(config::sample_cardinality);

		#pragma omp for reduction(+:reject)
		for (unsigned long j=0; j < config::runs; j++) {

			std::fill(freq_montecarlo.begin(), freq_montecarlo.end(), 0);

#if DEBUG
			if(omp_get_thread_num() == 0) {
				if ( j % (config::runs/80) == 0)
					std::cout << " " << j*8 << "/" << config::runs << std::endl;
			}
#endif
			MONTE_CARLO_GENERATOR(random_gen, freq_montecarlo);

			calculate_cumulative(freq_montecarlo, F_montecarlo);

			for (unsigned int i=0; i < config::sample_cardinality; i++) {
#if DEBUG
				std::cout << i << ": " << F_evt[i] << "==" <<  F_montecarlo[i] << std::endl;
#endif

				double diff = std::abs(F_evt[i] - F_montecarlo[i]);
				if (diff > ks_critical_value) {
					reject++;
					break;
				}
			}


		}
	
	}


	return reject;

}

template<int alpha_num, int alpha_den>
int run() noexcept {
	constexpr double alpha = ((double)alpha_num)/((double)alpha_den);

	#if TEST == TEST_KS
		unsigned long reject=perform_ks_run<alpha_num, alpha_den>();
	#elif TEST == TEST_AD
		unsigned long reject=perform_ad_run<alpha_num, alpha_den>();
	#else
		#error "Unknown test"
	#endif	

	unsigned long not_reject=0;
	not_reject = config::runs - reject;

	std::cout << std::setprecision(2) << alpha << "," << config::sample_cardinality << "," << reject << "," << not_reject << "," << std::setprecision(51) << ((double)reject) / (not_reject+reject) << std::endl;

	return not_reject;

}

int main(int argc, char* argv[]) {

	std::vector<unsigned int> cardinalities({50, 100, 250, 500, 750, 1000, 2500, 5000, 10000});

	for (const auto i : cardinalities) {
		config::sample_cardinality = i;
		int not_reject = run<5,100>();
		if (not_reject == 0) {
			break;
		}
	}

	for (const auto i : cardinalities) {
		config::sample_cardinality = i;
		int not_reject = run<1,100>();
		if (not_reject == 0) {
			break;
		}
	}

	return 0;
}




