if(Qt5_FOUND AND WIN32)
	add_executable(splittingVisualiser splittingVisualiser.cpp splittingVisualiser.h main.cpp)
	target_link_libraries(splittingVisualiser PUBLIC Qt5::Core Qt5::Gui Qt5::Widgets)
	target_link_libraries(splittingVisualiser PUBLIC networkReliabilityCommon mpfr ${Boost_LIBRARIES})
	target_compile_features(splittingVisualiser PUBLIC cxx_rvalue_references cxx_lambdas)
	target_include_directories(splittingVisualiser PUBLIC ${Boost_INCLUDE_DIRS})
endif()

