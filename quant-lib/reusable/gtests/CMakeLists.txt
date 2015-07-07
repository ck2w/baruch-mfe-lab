cmake_minimum_required(VERSION 2.6)
 
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ../ ../../utility/)
  
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests ../time_dependent_parameters.cpp ./time_dependent_parameters.t.cpp ../../utility/romberg_integrator.cpp) 
target_link_libraries(runTests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARY} pthread)
