cmake_minimum_required(VERSION 2.6)
 
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ../ ../../reusable/ ../../utility/)
  
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests ../bond.cpp ../../reusable/parameters.cpp ../../reusable/constant_parameters.cpp ../../reusable/time_dependent_parameters.cpp ../../utility/romberg_integrator.cpp ./bond.t.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARY} pthread)
