#ifndef DISTRIBUTIONS_HPP_
#define DISTRIBUTIONS_HPP_

#include <cmath>
#include <vector>

extern void fill_evt_cumulative(double xi, double mu, double sigma, std::vector<double> &output);
extern void fill_evt_cumulative_xi0(double mu, double sigma, std::vector<double> &output);

#endif
