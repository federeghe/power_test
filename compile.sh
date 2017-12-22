#!/bin/bash


build() {
	cd build
	make clean
	cmake -D CMAKE_CXX_FLAGS="-DH0=$1 -DH1=$2" ..
	make -j4
	mv power_test ../../power_test_$3
	cd ..
}

build 4 1 "h0_evt0.5_h1_t_student"
build 4 2 "h0_evt0.5_h1_uniform"
build 4 3 "h0_evt0.5_h1_evt0"
build 4 5 "h0_evt0.5_h1_evt-0.5"
build 4 6 "h0_evt0.5_h1_normal"

build 5 1 "h0_evt-0.5_h1_t_student"
build 5 2 "h0_evt-0.5_h1_uniform"
build 5 3 "h0_evt-0.5_h1_evt0"
build 5 5 "h0_evt-0.5_h1_evt-0.5"
build 5 6 "h0_evt-0.5_h1_normal"
