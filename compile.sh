#options="-static"
 options=" "

#optim="-Ofast"
 optim="-O3"

g++ src/*.cpp -c -std=c++11 $options $optim $memor

g++ -o demo *.o -std=c++11  $options $optim $memor

rm *.o
