#! /bin/bash

g++ -std=c++17 -O0 -DALG1 test_normalize_angle.cc -o nang10.out
g++ -std=c++17 -O0 -DALG2 test_normalize_angle.cc -o nang20.out
g++ -std=c++17 -O0 -DALG3 test_normalize_angle.cc -o nang30.out
g++ -std=c++17 -O0 -DALG4 test_normalize_angle.cc -o nang40.out

g++ -std=c++17 -O2 -DALG1 test_normalize_angle.cc -o nang12.out
g++ -std=c++17 -O2 -DALG2 test_normalize_angle.cc -o nang22.out
g++ -std=c++17 -O2 -DALG3 test_normalize_angle.cc -o nang32.out
g++ -std=c++17 -O2 -DALG4 test_normalize_angle.cc -o nang42.out
