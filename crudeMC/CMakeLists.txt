cmake_minimum_required(VERSION 3.1)
add_executable(crudeMC crudeMC.cpp)
target_link_libraries(crudeMC PUBLIC networkReliabilityCommon mpfr)
target_link_libraries(crudeMC PUBLIC ${Boost_LIBRARIES})
