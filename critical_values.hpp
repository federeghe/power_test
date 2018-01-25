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

template <int alpha_num, int alpha_den, bool upper_only=false>
class AndersonDarlingCV_UPPER {

	// From: "Approximation of modified Anderson-Darling test statisitcs for extreme value distributions with unknown shape parameter"
	//		J. Heo, H. Shin, W. Nam, J. Om, C. Jeong

};

template <>
class AndersonDarlingCV_UPPER<5, 100> {
public:
	static double get_critical_value(unsigned int cardinality, double shape_param) {
		switch(cardinality) {
		case 50:
			 if (shape_param==-0.5) return 1.30227776;
			 if (shape_param== 0.0) return 1.30199716;
			 if (shape_param== 0.5) return 1.306991668;
		break;
		case 100:
			 if (shape_param==-0.5) return 1.301841998;
			 if (shape_param== 0.0) return 1.306409406;
			 if (shape_param== 0.5) return 1.304260943;
		break;
		case 150:
			 if (shape_param==-0.5) return 1.303152793;
			 if (shape_param== 0.0) return 1.305088182;
			 if (shape_param== 0.5) return 1.301583148;
		break;
		case 200:
			 if (shape_param==-0.5) return 1.302306116;
			 if (shape_param== 0.0) return 1.303087216;
			 if (shape_param== 0.5) return 1.303143223;
		break;
		case 250:
			 if (shape_param==-0.5) return 1.299707246;
			 if (shape_param== 0.0) return 1.3007512;
			 if (shape_param== 0.5) return 1.302883938;
		break;
		case 300:
			 if (shape_param==-0.5) return 1.303490592;
			 if (shape_param== 0.0) return 1.302240504;
			 if (shape_param== 0.5) return 1.305810346;
		break;
		case 400:
			 if (shape_param==-0.5) return 1.30488691;
			 if (shape_param== 0.0) return 1.305784618;
			 if (shape_param== 0.5) return 1.303298465;
		break;
		case 500:
			 if (shape_param==-0.5) return 1.297476928;
			 if (shape_param== 0.0) return 1.301923002;
			 if (shape_param== 0.5) return 1.300142947;
		break;
		case 750:
			 if (shape_param==-0.5) return 1.303292006;
			 if (shape_param== 0.0) return 1.306332535;
			 if (shape_param== 0.5) return 1.298857487;
		break;
		case 1000:
			 if (shape_param==-0.5) return 1.299657242;
			 if (shape_param== 0.0) return 1.305620808;
			 if (shape_param== 0.5) return 1.305138494;
		break;
		case 2500:
			 if (shape_param==-0.5) return 1.304225486;
			 if (shape_param== 0.0) return 1.299562144;
			 if (shape_param== 0.5) return 1.302144513;
		break;
		case 5000:
			 if (shape_param==-0.5) return 1.301739083;
			 if (shape_param== 0.0) return 1.303781991;
			 if (shape_param== 0.5) return 1.304817548;
		break;
		case 10000:
			 if (shape_param==-0.5) return 1.302808865;
			 if (shape_param== 0.0) return 1.300236503;
			 if (shape_param== 0.5) return 1.299472353;
		break;
		}
		assert(false);
	}
};
template <>
class AndersonDarlingCV_UPPER<1, 100> {
public:
	static double get_critical_value(unsigned int cardinality, double shape_param) {
		switch(cardinality) {
		case 50:
			 if (shape_param==-0.5) return 2.079868153;
			 if (shape_param== 0.0) return 2.079540364;
			 if (shape_param== 0.5) return 2.080597454;
		break;
		case 100:
			 if (shape_param==-0.5) return 2.072568519;
			 if (shape_param== 0.0) return 2.070659818;
			 if (shape_param== 0.5) return 2.073940917;
		break;
		case 150:
			 if (shape_param==-0.5) return 2.073741786;
			 if (shape_param== 0.0) return 2.070340229;
			 if (shape_param== 0.5) return 2.06721026;
		break;
		case 200:
			 if (shape_param==-0.5) return 2.071085263;
			 if (shape_param== 0.0) return 2.061579587;
			 if (shape_param== 0.5) return 2.057173657;
		break;
		case 250:
			 if (shape_param==-0.5) return 2.071074206;
			 if (shape_param== 0.0) return 2.059420721;
			 if (shape_param== 0.5) return 2.056609281;
		break;
		case 300:
			 if (shape_param==-0.5) return 2.06048051;
			 if (shape_param== 0.0) return 2.065203622;
			 if (shape_param== 0.5) return 2.063771789;
		break;
		case 400:
			 if (shape_param==-0.5) return 2.060209887;
			 if (shape_param== 0.0) return 2.066727898;
			 if (shape_param== 0.5) return 2.060353413;
		break;
		case 500:
			 if (shape_param==-0.5) return 2.055336251;
			 if (shape_param== 0.0) return 2.057690045;
			 if (shape_param== 0.5) return 2.053739255;
		break;
		case 750:
			 if (shape_param==-0.5) return 2.0591838;
			 if (shape_param== 0.0) return 2.059023775;
			 if (shape_param== 0.5) return 2.065388522;
		break;
		case 1000:
			 if (shape_param==-0.5) return 2.057976117;
			 if (shape_param== 0.0) return 2.057185359;
			 if (shape_param== 0.5) return 2.063604049;
		break;
		case 2500:
			 if (shape_param==-0.5) return 2.051217776;
			 if (shape_param== 0.0) return 2.059171445;
			 if (shape_param== 0.5) return 2.061884222;
		break;
		case 5000:
			 if (shape_param==-0.5) return 2.060789138;
			 if (shape_param== 0.0) return 2.052528229;
			 if (shape_param== 0.5) return 2.059770473;
		break;
		case 10000:
			 if (shape_param==-0.5) return 2.054734726;
			 if (shape_param== 0.0) return 2.062242055;
			 if (shape_param== 0.5) return 2.064245176;
		break;
		}
		assert(false);
	}
};
