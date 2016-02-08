cmake_minimum_required(VERSION 2.6)
 
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ../ ../../reusable/ ../../utility/)
  
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests ../european_option.cpp ../payoff.cpp ../vanilla_option.cpp ../european_payoff.cpp ../../reusable/constant_parameters.cpp ../../utility/cumulative_normal.cpp ./black_scholes.t.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARY} pthread)
