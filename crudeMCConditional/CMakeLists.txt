cmake_minimum_required(VERSION 3.1)
add_executable(crudeMCConditional crudeMCConditional.cpp)
target_link_libraries(crudeMCConditional PUBLIC networkReliabilityCommon mpfr ${Boost_LIBRARIES})
