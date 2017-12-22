#include "critical_values.hpp"
#include "montecarlo.hpp"
#include "distributions.hpp"
#include "simulation.hpp"

#include <iomanip>    
#include <iostream>
#include <omp.h>


#define DEBUG 0


unsigned int config::sample_cardinality;

void run(double xi) noexcept {

	std::vector<double> F_evt_0;
	std::vector<double> F_evt_eps;

	F_evt_0.resize(config::size);
	F_evt_eps.resize(config::size);

	fill_evt_cumulative_xi0(0, 1, F_evt_0);
	fill_evt_cumulative(xi, 0, 1, F_evt_eps);

	double ks_critical_value1 = get_ks_critical<1, 100>(config::sample_cardinality);
	double ks_critical_value2 = get_ks_critical<5, 100>(config::sample_cardinality);

#if DEBUG
	std::cout << i << ": " << F_evt_0[i] << "==" <<  F_evt_eps[i] << std::endl;
#endif

	std::cout << config::sample_cardinality << " " << xi << " ";

	bool reject1 = false;
	bool reject2 = false;

	for (unsigned int i=0; i < config::size; i++) {
		double diff = std::abs(F_evt_0[i] - F_evt_eps[i]);
		if (diff > ks_critical_value1) {
			reject1=true;
		}
		if (diff > ks_critical_value2) {
			reject2=true;
		}

	}

	if (reject1) {
		std::cout << "1";	// Reject
	} 
	else {
		std::cout << "0";	// Not reject
	}

	std::cout << " ";

	if (reject2) {
		std::cout << "1";	// Reject
	} 
	else {
		std::cout << "0";	// Not reject
	}


	std::cout << std::endl;

}

int main(int argc, char* argv[]) {



	std::vector<unsigned int> cardinalities;
	for (int i=0; i<1000; i++) {
		cardinalities.push_back((i+1)*10);
	}

	for (double xi=0.001; xi <= 0.5; xi += 0.001) {
		for (const auto i : cardinalities) {
			config::sample_cardinality = i;
			run(xi);
		}

	}
	
	return 0;
}




