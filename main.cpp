#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"

#include <iomanip>    
#include <iostream>
#include <fstream>
#include <omp.h>

unsigned int config::sample_cardinality;

std::ofstream out_file("h0_e0_h1_t_student.data");


template<int alpha_num, int alpha_den>
int run() noexcept {

	constexpr double alpha = ((double)alpha_num)/((double)alpha_den);

	std::vector<double> F_evt;

	F_evt.resize(config::size);

	fill_evt_cumulative_xi0(0, 1, F_evt);


	double ks_critical_value = get_kd_critical<alpha_num, alpha_den>(config::sample_cardinality);

	unsigned long reject=0;
	
	#pragma omp parallel firstprivate(F_evt)
	{

		thread_local static std::mt19937 random_gen(std::random_device{}());


		std::vector<unsigned int> freq_montecarlo;
		std::vector<double> F_montecarlo;
		freq_montecarlo.resize(config::size);
		F_montecarlo.resize(config::size);

		#pragma omp for reduction(+:reject)
		for (unsigned long j=0; j < config::runs; j++) {

			std::fill(freq_montecarlo.begin(), freq_montecarlo.end(), 0);

/*			if(omp_get_thread_num() == 0) {
				if ( j % (config::runs/80) == 0)
					std::cout << " " << j*8 << "/" << config::runs << std::endl;
			}
*/
//			montecarlo_frequencies<T_STUDENT>(random_gen, freq_montecarlo);

			montecarlo_frequencies_evt<0,1>(random_gen, freq_montecarlo);

			calculate_cumulative(freq_montecarlo, F_montecarlo);

			for (unsigned int i=0; i < config::size; i++) {

				double diff = std::abs(F_evt[i] - F_montecarlo[i]);
				if (diff > ks_critical_value) {
					reject++;
					break;
				}
			}

		}
	
	}

	unsigned long not_reject=0;
	not_reject = config::runs - reject;

	out_file << alpha << "," << config::sample_cardinality << "," << reject << "," << not_reject << "," << std::setprecision(51) << ((double)reject) / (not_reject+reject) << std::endl;


	std::cout << "REJECT: " << reject << " NOT REJECT: " << not_reject << " POWER: " << std::setprecision(51) << ((double)reject) / (not_reject+reject) << std::endl;

	return not_reject;

}

int main(int argc, char* argv[]) {



	std::vector<unsigned int> cardinalities({50, 100, 250, 500, 750, 1000, 2500, 5000, 10000});

	std::cout << std::endl << "############################################" << std::endl;
	std::cout << "####  ALPHA=0.05" << std::endl;
	std::cout << "############################################" << std::endl << std::endl;

	for (const auto i : cardinalities) {
		std::cout << "-------------------------" << std::endl;
		std::cout << "N=" << i << std::endl;
		config::sample_cardinality = i;
		int not_reject = run<5,100>();
		if (not_reject == 0) {
			break;
		}

	}

	std::cout << std::endl << "############################################" << std::endl;
	std::cout << "####  ALPHA=0.01" << std::endl;
	std::cout << "############################################" << std::endl << std::endl;


	for (const auto i : cardinalities) {
		std::cout << "-------------------------" << std::endl;
		std::cout << "N=" << i << std::endl;
		config::sample_cardinality = i;
		int not_reject = run<1,100>();
		if (not_reject == 0) {
			break;
		}

	}


	return 0;
}




