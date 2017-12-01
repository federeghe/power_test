cmake_minimum_required(VERSION 2.8)
project(power_test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -lpthread -Wall -pedantic")


add_executable(power_test critical_values.cpp distributions.cpp montecarlo.cpp main.cpp)

