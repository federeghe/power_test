#include <cassert>
#include <cmath>


template <int alpha_num, int alpha_den>
inline double get_ks_critical(unsigned long n) {
	constexpr double alpha = ((double)alpha_num)/((double)alpha_den);
	return std::sqrt(-0.5 * std::log(alpha/2))/std::sqrt(n);

	// Accuracy good for n >= 30
	

}



template <int alpha_num, int alpha_den>
class AndersonDarlingCV_EV {

	// From: "Approximation of modified Anderson-Darling test statisitcs for extreme value distributions with unknown shape parameter"
	//		J. Heo, H. Shin, W. Nam, J. Om, C. Jeong

};

template <>
class AndersonDarlingCV_EV<5, 100> {

public:

	static constexpr double safety_margin = 0.05;
	
	static double get_critical_value(unsigned long n, double shape_param) noexcept {

		constexpr double a = 0.5507;
		constexpr double b = 0.7928;
		constexpr double c = -3.5521;
		constexpr double d = -0.2103;
		constexpr double e = 0.3388;
		

		double table_cv = a + b / n + c / (n*n) + d * shape_param + e * shape_param * shape_param;

		return table_cv + table_cv * safety_margin;
	}

};

template <>
class AndersonDarlingCV_EV<1, 100> {

public:
	static constexpr double safety_margin = 0.05;
	
	static double get_critical_value(unsigned long n, double shape_param) noexcept {

		constexpr double a = 0.7385;
		constexpr double b = 2.5585;
		constexpr double c = -14.2247;
		constexpr double d = -0.5708;
		constexpr double e = 0.8825;
		

		double table_cv = a + b / n + c / (n*n) + d * shape_param + e * shape_param * shape_param;

		return table_cv + table_cv * safety_margin;
	}

};
