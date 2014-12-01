cmake_minimum_required(VERSION 3.1)
add_executable(exhaustiveProbability main.cpp)
target_link_libraries(exhaustiveProbability PUBLIC networkReliabilityCommon mpfr)
target_link_libraries(exhaustiveProbability PUBLIC ${Boost_LIBRARIES})
