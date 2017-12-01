#include <cmath>


template <int alpha_num, int alpha_den>
inline double get_kd_critical(unsigned long n) {
	constexpr double alpha = ((double)alpha_num)/((double)alpha_den);
	return std::sqrt(-0.5 * std::log(alpha/2))/std::sqrt(n);
}

