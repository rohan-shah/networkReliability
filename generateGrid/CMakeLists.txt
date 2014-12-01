cmake_minimum_required(VERSION 3.1)
add_executable(generateGrid generateGrid.cpp)
target_link_libraries(generateGrid PUBLIC networkReliabilityCommon mpfr ${Boost_LIBRARIES})
