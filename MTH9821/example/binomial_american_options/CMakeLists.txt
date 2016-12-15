cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-long-long -Wno-deprecated -pedantic")

include_directories(./ ../../binomial_tree/ ../../options/)
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runBinomialAmerican
               ./main.cpp
               ../../binomial_tree/binomial.cpp
               ../../options/black_scholes.cpp)

