#include <cassert>
#include <cmath>

template <int alpha_num, int alpha_den>
inline double get_ks_critical(unsigned long n) {
	constexpr double alpha = ((double)alpha_num)/((double)alpha_den);
	return std::sqrt(-0.5 * std::log(alpha/2))/std::sqrt(n);

	// Accuracy good for n >= 30
	

}



template <int alpha_num, int alpha_den, bool upper_only=false>
class AndersonDarlingCV {

	// From: "Approximation of modified Anderson-Darling test statisitcs for extreme value distributions with unknown shape parameter"
	//		J. Heo, H. Shin, W. Nam, J. Om, C. Jeong

};

template <>
class AndersonDarlingCV<5, 100> {
public:
	static double get_critical_value(unsigned int cardinality, double shape_param) {
		switch(cardinality) {
		case 50:
			 if (shape_param==-0.5) return 2.47412;
			 if (shape_param== 0.0) return 2.50976;
			 if (shape_param== 0.5) return 2.53496;
		break;
		case 100:
			 if (shape_param==-0.5) return 2.49509;
			 if (shape_param== 0.0) return 2.52714;
			 if (shape_param== 0.5) return 2.42694;
		break;
		case 250:
			 if (shape_param==-0.5) return 2.48027;
			 if (shape_param== 0.0) return 2.52288;
			 if (shape_param== 0.5) return 2.50464;
		break;
		case 500:
			 if (shape_param==-0.5) return 2.46681;
			 if (shape_param== 0.0) return 2.44137;
			 if (shape_param== 0.5) return 2.48736;
		break;
		case 750:
			 if (shape_param==-0.5) return 2.50228;
			 if (shape_param== 0.0) return 2.48115;
			 if (shape_param== 0.5) return 2.50141;
		break;
		case 1000:
			 if (shape_param==-0.5) return 2.44515;
			 if (shape_param== 0.0) return 2.46962;
			 if (shape_param== 0.5) return 2.52046;
		break;
		case 2500:
			 if (shape_param==-0.5) return 2.56467;
			 if (shape_param== 0.0) return 2.54922;
			 if (shape_param== 0.5) return 2.44433;
		break;
		case 5000:
			 if (shape_param==-0.5) return 2.47965;
			 if (shape_param== 0.0) return 2.48939;
			 if (shape_param== 0.5) return 2.47019;
		break;
		case 10000:
			 if (shape_param==-0.5) return 2.51519;
			 if (shape_param== 0.0) return 2.53523;
			 if (shape_param== 0.5) return 2.5269;
		break;
		}
		assert(false);
	}
};
template <>
class AndersonDarlingCV<1, 100> {
public:
	static double get_critical_value(unsigned int cardinality, double shape_param) {
		switch(cardinality) {
		case 50:
			 if (shape_param==-0.5) return 3.77154;
			 if (shape_param== 0.0) return 3.86773;
			 if (shape_param== 0.5) return 3.74202;
		break;
		case 100:
			 if (shape_param==-0.5) return 3.86384;
			 if (shape_param== 0.0) return 3.87911;
			 if (shape_param== 0.5) return 3.89836;
		break;
		case 250:
			 if (shape_param==-0.5) return 3.81976;
			 if (shape_param== 0.0) return 3.93195;
			 if (shape_param== 0.5) return 3.93606;
		break;
		case 500:
			 if (shape_param==-0.5) return 3.92627;
			 if (shape_param== 0.0) return 3.86792;
			 if (shape_param== 0.5) return 3.84262;
		break;
		case 750:
			 if (shape_param==-0.5) return 3.90418;
			 if (shape_param== 0.0) return 3.68703;
			 if (shape_param== 0.5) return 3.86732;
		break;
		case 1000:
			 if (shape_param==-0.5) return 3.91767;
			 if (shape_param== 0.0) return 3.88882;
			 if (shape_param== 0.5) return 3.8788;
		break;
		case 2500:
			 if (shape_param==-0.5) return 4.06825;
			 if (shape_param== 0.0) return 3.98086;
			 if (shape_param== 0.5) return 3.92626;
		break;
		case 5000:
			 if (shape_param==-0.5) return 3.8144;
			 if (shape_param== 0.0) return 3.88781;
			 if (shape_param== 0.5) return 3.89555;
		break;
		case 10000:
			 if (shape_param==-0.5) return 3.74082;
			 if (shape_param== 0.0) return 3.80453;
			 if (shape_param== 0.5) return 3.88444;
		break;
		}
		assert(false);
	}
};

