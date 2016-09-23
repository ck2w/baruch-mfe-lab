cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-long-long -Wno-deprecated -pedantic")

include_directories(./
                    ../load/
                    ../cholesky 
                    ../triangular_solve
                    ../band_op
                    ../../../lib/eigen/)
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(run
               ./main.cpp
               ./linear_solve.cpp
               ../load/load.cpp
               ../cholesky/cholesky.cpp
               ../band_op/band_conversion.cpp
               ../triangular_solve/triangular_solve.cpp)

