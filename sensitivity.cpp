#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"

#include <iomanip>
#include <iostream>
#include <omp.h>

#define NR_RUNS 1000
#define DEBUG 0


unsigned int config::sample_cardinality;


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
	for (int i=0; i<1000; i++) {
		cardinalities.push_back((i+1)*10);
	}


	for (const auto i : cardinalities) {
		config::sample_cardinality = i;
		static_for<1,501>()();
	}
	
	return 0;
}




