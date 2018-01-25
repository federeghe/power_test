#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_


#define DIST_T_STUDENT 1
#define DIST_UNIFORM 2
#define DIST_EVT0 3
#define DIST_EVT_0_5 4
#define DIST_EVT_N_0_5 5
#define DIST_NORMAL 6

//#define H1 DIST_EVT_0_5

//#define H0 DIST_EVT_0_5

//#define DEBUG 1

#define TEST_KS 1
#define TEST_AD 2
#define TEST_MAD 3

#if H0 == DIST_EVT0
	#define EVT_PARAM_NUM 0
	#define EVT_PARAM_DEN 1
constexpr double evt_param = 0.0;
#elif H0 == DIST_EVT_0_5
	#define EVT_PARAM_NUM 1
	#define EVT_PARAM_DEN 2
constexpr double evt_param = 0.5;
#elif H0 == DIST_EVT_N_0_5
	#define EVT_PARAM_NUM -1
	#define EVT_PARAM_DEN 2
constexpr double evt_param = -0.5;
#else
	#error "Unknown distribution"
#endif

#if H1 == DIST_T_STUDENT
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<T_STUDENT>
	#define MONTE_CARLO_SAMPLING montecarlo_t_student_sample
#elif H1 == DIST_UNIFORM
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<UNIFORM>
	#define MONTE_CARLO_SAMPLING montecarlo_uniform_sample
#elif H1 == DIST_NORMAL
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies<NORMAL>
	#define MONTE_CARLO_SAMPLING montecarlo_normal_sample
#elif H1 == DIST_EVT0
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<0,1>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<0, 1>
#elif H1 == DIST_EVT_0_5
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<1,2>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<1, 2>
#elif H1 == DIST_EVT_N_0_5
	#define MONTE_CARLO_GENERATOR montecarlo_frequencies_evt<-1,2>
	#define MONTE_CARLO_SAMPLING montecarlo_evt_sample<-1, 2>
#else
	#error "Unknown"
#endif


namespace config {
	constexpr double dist_min=-15;
	constexpr double dist_max=15;
	constexpr unsigned long runs = 1e9;
	extern unsigned int sample_cardinality;
}

#endif
