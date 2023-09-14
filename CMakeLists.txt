cmake_minimum_required(VERSION 3.1)

project( FastContainerTest) 
set(exec_name ${PROJECT_NAME}.bin)

include_directories( "include" )

set( CMAKE_CXX_STANDARD 17 )
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

aux_source_directory( ./src sources)
add_executable( ${exec_name} test.cpp ${sources})
target_link_libraries( ${exec_name} )

install( TARGETS ${exec_name} DESTINATION ${CMAKE_SOURCE_DIR}/bin )
