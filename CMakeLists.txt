cmake_minimum_required(VERSION 3.1)

project( FastContainerTest) 
set(exec_name ${PROJECT_NAME}.bin)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT 6.0)

include_directories( "include" )
include_directories( ${ROOT_INCLUDE_DIRS} )

set( CMAKE_CXX_STANDARD 17 )
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

aux_source_directory( ./src sources)
add_executable( ${exec_name} test.cpp ${sources})
target_link_libraries( ${exec_name} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES})

install( TARGETS ${exec_name} DESTINATION ${CMAKE_SOURCE_DIR}/bin )
