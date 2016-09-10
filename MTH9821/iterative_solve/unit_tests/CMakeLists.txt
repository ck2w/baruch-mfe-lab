cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-long-long -Wno-deprecated -pedantic")

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} 
                    ../ 
                    ../../../../lib/eigen/
                    ../../triangular_solve/)
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests
               ./iterative_solve.t.cpp
               ../linear_iterate.cpp
               ../sor.cpp
               ../jacobi.cpp
               ../../triangular_solve/triangular_solve.cpp)

target_link_libraries(runTests 
                      ${GTEST_LIBRARIES} 
                      ${GTEST_MAIN_LIBRARY} 
                      pthread)
