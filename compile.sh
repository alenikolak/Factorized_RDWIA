#bin/bash

#optim="-Ofast"
optim="-O3"

g++ src/Crosssections.cpp src/Formfactors.cpp src/Hadronic.cpp src/Leptonic.cpp src/S_tables.cpp -c -std=c++11 $optim

g++ -o benchmarks src/Benchmarks.cpp *.o -std=c++11  $optim
g++ -o demo src/Demo.cpp *.o -std=c++11  $optim

rm *.o
