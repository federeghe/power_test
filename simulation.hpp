#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

namespace config {
	constexpr double dist_min=-15;
	constexpr double dist_max=15;
	constexpr double step=0.01;
	constexpr unsigned int size = (dist_max - dist_min) / step;
	constexpr unsigned long runs = 1e7;
	extern unsigned int sample_cardinality;
}

#endif
