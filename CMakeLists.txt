cmake_minimum_required(VERSION 2.8)
project(power_test)

SET(CMAKE_C_COMPILER mpiicc)
SET(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -lpthread -Wall -pedantic -std=c++11")


add_executable(power_test critical_values.cpp distributions.cpp montecarlo.cpp main.cpp)
add_executable(sensitivity_test distributions.cpp montecarlo.cpp sensitivity.cpp)
add_executable(ad_critical_value ad_critical_value.cpp distributions.cpp montecarlo.cpp)

target_link_libraries(power_test mpi)


