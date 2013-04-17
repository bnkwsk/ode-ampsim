#!/bin/bash
maxima < algebra/algebra.max
g++ -I/home/prezes/lib/flens/ -std=c++11 -O3 symbolic_main.cpp -o inz -lsndfile -lcln
