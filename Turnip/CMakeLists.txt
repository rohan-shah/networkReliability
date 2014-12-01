cmake_minimum_required(VERSION 3.1)
add_executable(turnip turnip.cpp)
target_link_libraries(turnip PUBLIC networkReliabilityCommon mpfr ${Boost_LIBRARIES})
