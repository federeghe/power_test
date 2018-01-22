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
			 if (shape_param==-0.5) return 2.491652034;
			 if (shape_param== 0.0) return 2.493387134;
			 if (shape_param== 0.5) return 2.499743513;
		break;
		case 100:
			 if (shape_param==-0.5) return 2.491294015;
			 if (shape_param== 0.0) return 2.495910411;
			 if (shape_param== 0.5) return 2.498130845;
		break;
		case 150:
			 if (shape_param==-0.5) return 2.492142674;
			 if (shape_param== 0.0) return 2.497623943;
			 if (shape_param== 0.5) return 2.495339529;
		break;
		case 200:
			 if (shape_param==-0.5) return 2.491183585;
			 if (shape_param== 0.0) return 2.492421764;
			 if (shape_param== 0.5) return 2.494171316;
		break;
		case 250:
			 if (shape_param==-0.5) return 2.494747063;
			 if (shape_param== 0.0) return 2.490413056;
			 if (shape_param== 0.5) return 2.491140562;
		break;
		case 300:
			 if (shape_param==-0.5) return 2.496596727;
			 if (shape_param== 0.0) return 2.496531505;
			 if (shape_param== 0.5) return 2.493510858;
		break;
		case 400:
			 if (shape_param==-0.5) return 2.490163161;
			 if (shape_param== 0.0) return 2.494845986;
			 if (shape_param== 0.5) return 2.497766782;
		break;
		case 500:
			 if (shape_param==-0.5) return 2.494244971;
			 if (shape_param== 0.0) return 2.49511289;
			 if (shape_param== 0.5) return 2.494650458;
		break;
		case 750:
			 if (shape_param==-0.5) return 2.49998025;
			 if (shape_param== 0.0) return 2.495524615;
			 if (shape_param== 0.5) return 2.494604129;
		break;
		case 1000:
			 if (shape_param==-0.5) return 2.48901915;
			 if (shape_param== 0.0) return 2.487947993;
			 if (shape_param== 0.5) return 2.489190879;
		break;
		case 2500:
			 if (shape_param==-0.5) return 2.492478225;
			 if (shape_param== 0.0) return 2.494547674;
			 if (shape_param== 0.5) return 2.498825044;
		break;
		case 5000:
			 if (shape_param==-0.5) return 2.493096296;
			 if (shape_param== 0.0) return 2.49506778;
			 if (shape_param== 0.5) return 2.495002774;
		break;
		case 10000:
			 if (shape_param==-0.5) return 2.49291653;
			 if (shape_param== 0.0) return 2.480720271;
			 if (shape_param== 0.5) return 2.494940591;
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
			 if (shape_param==-0.5) return 3.881350134;
			 if (shape_param== 0.0) return 3.891467406;
			 if (shape_param== 0.5) return 3.892237753;
		break;
		case 100:
			 if (shape_param==-0.5) return 3.878370792;
			 if (shape_param== 0.0) return 3.886309477;
			 if (shape_param== 0.5) return 3.874740472;
		break;
		case 150:
			 if (shape_param==-0.5) return 3.884910892;
			 if (shape_param== 0.0) return 3.884027279;
			 if (shape_param== 0.5) return 3.881346088;
		break;
		case 200:
			 if (shape_param==-0.5) return 3.876058408;
			 if (shape_param== 0.0) return 3.875536532;
			 if (shape_param== 0.5) return 3.891019076;
		break;
		case 250:
			 if (shape_param==-0.5) return 3.872636202;
			 if (shape_param== 0.0) return 3.891804911;
			 if (shape_param== 0.5) return 3.8941785;
		break;
		case 300:
			 if (shape_param==-0.5) return 3.854080842;
			 if (shape_param== 0.0) return 3.886549965;
			 if (shape_param== 0.5) return 3.870747601;
		break;
		case 400:
			 if (shape_param==-0.5) return 3.889208171;
			 if (shape_param== 0.0) return 3.884914412;
			 if (shape_param== 0.5) return 3.88164018;
		break;
		case 500:
			 if (shape_param==-0.5) return 3.886266838;
			 if (shape_param== 0.0) return 3.867722582;
			 if (shape_param== 0.5) return 3.885585846;
		break;
		case 750:
			 if (shape_param==-0.5) return 3.870929838;
			 if (shape_param== 0.0) return 3.883727267;
			 if (shape_param== 0.5) return 3.868245583;
		break;
		case 1000:
			 if (shape_param==-0.5) return 3.882895784;
			 if (shape_param== 0.0) return 3.882021218;
			 if (shape_param== 0.5) return 3.889765691;
		break;
		case 2500:
			 if (shape_param==-0.5) return 3.886627592;
			 if (shape_param== 0.0) return 3.878939207;
			 if (shape_param== 0.5) return 3.890788073;
		break;
		case 5000:
			 if (shape_param==-0.5) return 3.881622021;
			 if (shape_param== 0.0) return 3.893294697;
			 if (shape_param== 0.5) return 3.868217466;
		break;
		case 10000:
			 if (shape_param==-0.5) return 3.886011819;
			 if (shape_param== 0.0) return 3.855781773;
			 if (shape_param== 0.5) return 3.877462098;
		break;
		}
		assert(false);
	}
};

