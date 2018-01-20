#!/bin/bash


build() {
	cd build
	make clean
	cmake -D CMAKE_CXX_FLAGS="-DH0=$1 -DH1=$2 -DTEST=$3" ..
	make -j4
	mv power_test ../../power_test_${4}_statistic_${3}
	cd ..
}

build 3 1 2 "h0_evt0_h1_t_student"
build 3 2 2 "h0_evt0_h1_uniform"
build 3 4 2 "h0_evt0_h1_evt0.5"
build 3 5 2 "h0_evt0_h1_evt-0.5"
build 3 6 2 "h0_evt0_h1_normal"

build 4 1 2 "h0_evt0.5_h1_t_student"
build 4 2 2 "h0_evt0.5_h1_uniform"
build 4 3 2 "h0_evt0.5_h1_evt0"
build 4 5 2 "h0_evt0.5_h1_evt-0.5"
build 4 6 2 "h0_evt0.5_h1_normal"

build 5 1 2 "h0_evt-0.5_h1_t_student"
build 5 2 2 "h0_evt-0.5_h1_uniform"
build 5 3 2 "h0_evt-0.5_h1_evt0"
build 5 4 2 "h0_evt-0.5_h1_evt0.5"
build 5 6 2 "h0_evt-0.5_h1_normal"

