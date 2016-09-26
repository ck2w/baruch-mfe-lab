cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-long-long -Wno-deprecated -pedantic")

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ../)
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests
               ./binomial.t.cpp
               ../black_scholes.cpp
               ../binomial.cpp)

target_link_libraries(runTests 
                      ${GTEST_LIBRARIES} 
                      ${GTEST_MAIN_LIBRARY} 
                      pthread)
