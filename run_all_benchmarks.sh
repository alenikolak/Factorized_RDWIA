#bin/bash

echo "Reading the tables takes most of the time"
./benchmarks <<< "C12 1 1" > Benchmarks_C12.out
./benchmarks <<< "C12 1 2" >> Benchmarks_C12.out
./benchmarks <<< "C12 2 1" >> Benchmarks_C12.out
./benchmarks <<< "C12 2 2" >> Benchmarks_C12.out
echo "Half way there"
./benchmarks <<< "C12 4 1" >> Benchmarks_C12.out
./benchmarks <<< "C12 4 2" >> Benchmarks_C12.out
./benchmarks <<< "C12 5 1" >> Benchmarks_C12.out
./benchmarks <<< "C12 5 2" >> Benchmarks_C12.out
