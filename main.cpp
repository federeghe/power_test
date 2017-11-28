#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"

#include <iomanip>    
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[]) {

	rng.seed(std::random_device()());


	std::vector<double> F_evt;


	F_evt.resize(config::size);

	fill_evt_cumulative_xi0(0, 1, F_evt);


	double ks_critical_value = get_kd_critical(0.2, config::sample_cardinality);

	unsigned long reject=0;
	unsigned long not_reject=0;
	unsigned long counter=0;
	
	#pragma omp parallel firstprivate(F_evt)
	{

		std::vector<unsigned int> freq_montecarlo;
		std::vector<double> F_montecarlo;
		freq_montecarlo.resize(config::size);
		F_montecarlo.resize(config::size);

		#pragma omp for reduction(+:reject,not_reject)
		for (unsigned long j=0; j < config::runs; j++) {

			std::fill(freq_montecarlo.begin(), freq_montecarlo.end(), 0);

			if(omp_get_thread_num() == 0) {
				counter++;
				if ( counter % (config::runs/50) == 0)
					std::cout << (reject+not_reject) * omp_get_num_threads() << "/" << config::runs << std::endl;
			}

			montecarlo_frequencies(freq_montecarlo);

			calculate_cumulative(freq_montecarlo, F_montecarlo);

			for (unsigned int i=0; i<freq_montecarlo.size(); i++) {
			}

			double max_diff=0;
			for (unsigned int i=0; i<freq_montecarlo.size(); i++) {
				max_diff = std::max(max_diff, std::abs(F_evt[i] - F_montecarlo[i]));
			}

			if (max_diff > ks_critical_value) {
				reject++;
			} else {
				not_reject++;
			}
		}
	
	}

	std::cout << std::endl << "REJECT: " << reject << " NOT REJECT: " << not_reject << " POWER: " << std::setprecision(51) << ((double)reject) / (not_reject+reject) << std::endl;

	return 0;

}




