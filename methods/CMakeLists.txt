cmake_minimum_required(VERSION 3.1)
set(CppFiles crudeMCImpl.cpp)
set(HeaderFiles crudeMCImpl.h)

add_library(methods STATIC ${CppFiles} ${HeaderFiles})
target_link_libraries(methods PUBLIC networkReliabilityCommon) 

if(WIN32)
	set_target_properties(methods PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/x64/ PDB_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/x64/)
endif()
