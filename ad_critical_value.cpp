#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"
#include "test_statistics.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <omp.h>

#define RUNS 1000000

unsigned int config::sample_cardinality;

void compute(double s) {
	std::vector<double> crits_1;	// EVT<-0.5>
	std::vector<double> crits_2;	// EVT<0>
	std::vector<double> crits_3;	// EVT<0.5>

	crits_1.resize(RUNS);
	crits_2.resize(RUNS);
	crits_3.resize(RUNS);



	#pragma omp parallel 
	{
		std::vector<double> sample_1;
		std::vector<double> sample_2;
		std::vector<double> sample_3;
		sample_1.resize(config::sample_cardinality);
		sample_2.resize(config::sample_cardinality);
		sample_3.resize(config::sample_cardinality);

		std::mt19937 random_gen(std::random_device{}());
		#pragma omp for
		for (unsigned int i=0; i < RUNS; i++) {
			for (unsigned int j=0; j < config::sample_cardinality; j++) {
				sample_1[j] = montecarlo_evt_sample<-1, 2>(random_gen);
				sample_2[j] = montecarlo_evt_sample<0, 1> (random_gen);
				sample_3[j] = montecarlo_evt_sample<1, 2> (random_gen);
			}
			std::sort(sample_1.begin(), sample_1.end());
			std::sort(sample_2.begin(), sample_2.end());
			std::sort(sample_3.begin(), sample_3.end());
			crits_1[i] = get_ad_statistic<-1, 2>(sample_1);
			crits_2[i] = get_ad_statistic<0, 1> (sample_2);
			crits_3[i] = get_ad_statistic<1, 2> (sample_3);
		}

		#pragma omp sections
		{
			#pragma omp section
			std::sort(crits_1.begin(),crits_1.end());
			#pragma omp section
			std::sort(crits_2.begin(),crits_2.end());
			#pragma omp section
			std::sort(crits_3.begin(),crits_3.end());
		}
	}

	std::cout << "\t\tcase " << config::sample_cardinality << ":" << std::endl;
	std::cout << "\t\t\t if (shape_param==-0.5) return "	<< std::setprecision(10) << crits_1[(int)(s * RUNS)] << ';' << std::endl;
	std::cout << "\t\t\t if (shape_param== 0.0) return " 	<< std::setprecision(10) << crits_2[(int)(s * RUNS)] << ';' <<  std::endl;
	std::cout << "\t\t\t if (shape_param== 0.5) return "	<< std::setprecision(10) << crits_3[(int)(s * RUNS)] << ';' <<  std::endl;
	std::cout << "\t\tbreak;" << std::endl;

}

int main() {

	std::vector<unsigned int> cardinalities({50 /*, 100, 200, 250, 300, 400, 500, 750, 1000, 2500, 5000, 10000*/});

	std::vector<double> significance_level({0.95, 0.99});

	for (const auto s : significance_level) {
		std::cout << "template <>" << std::endl;
		std::cout << "class AndersonDarlingCV<" << (int)((1-s)*100) << ", 100> {" << std::endl;
		std::cout << "public:" << std::endl;

		std::cout << "\tstatic double get_critical_value(unsigned int cardinality, double shape_param) {" << std::endl;
		std::cout << "\t\tswitch(cardinality) {" << std::endl;

		for (const auto i : cardinalities) {
			config::sample_cardinality = i;
			compute(s);
		}
		std::cout << "\t\t}" << std::endl;
		std::cout << "\t\tassert(false);" << std::endl;
		std::cout << "\t}" << std::endl;
		std::cout << "};" << std::endl;
	}
	return 0;
}
