cmake_minimum_required(VERSION 2.6)

include_directories(./ ../load/ ../../../lib/eigen/)
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(run
               ./main.cpp
               ./triangular_solve.cpp
               ../load/load.cpp)

