cmake_minimum_required(VERSION 3.1)
add_executable(crudeMCConditionalFunction crudeMCConditionalFunction.cpp)

target_link_libraries(crudeMCConditionalFunction PUBLIC networkReliabilityCommon mpfr)
target_compile_definitions(crudeMCConditionalFunction PUBLIC BOOST_ALL_NO_LIB BOOST_ALL_DYN_LINK)
if(WIN32)
	set_target_properties(crudeMCConditionalFunction PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/x64/ PDB_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/x64/)
endif()
