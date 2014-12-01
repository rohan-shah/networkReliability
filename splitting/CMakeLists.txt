cmake_minimum_required(VERSION 3.1)
find_package(aliasMethod REQUIRED)
add_executable(splitting main.cpp)
target_link_libraries(splitting PUBLIC networkReliabilityCommon mpfr ${Boost_LIBRARIES} aliasMethod)
