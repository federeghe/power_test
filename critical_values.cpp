#include <cmath>

double get_kd_critical(double alpha, unsigned long n) {
	return std::sqrt(-0.5 * std::log(alpha/2))/std::sqrt(n);
}

