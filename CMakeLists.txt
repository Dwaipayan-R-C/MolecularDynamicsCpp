cmake_minimum_required(VERSION 3.11)
set( CMAKE_CXX_COMPILER "/usr/bin/g++" )
set( CMAKE_C_COMPILER "/usr/bin/gcc" )
project(yamd)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)
#Adding gtset
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#////////
add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(MYPROJECT_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/header/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/atom_structure.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/xyz.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/header/force.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/lj_direct_summation.h

        )

#        )

set(MYPROJECT_SRCS
        # Add source files here
#         ${CMAKE_CURRENT_SOURCE_DIR}/verlet.cpp
        )

set(MYPROJECT_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        header/xyz.h xyz.cpp header/verlet.h verlet.cpp header/lj_direct_summation.h lj_direct_summation.cpp
        lj_kinetic_potential.cpp header/lj_kinetic_and_potential.h header/berendsen_thermostat.h
        berendsen_thermostat.cpp header/neighbors.h neighbors.cpp header/gupta.h gupta.cpp)

#
#add_subdirectory(tests)


include(GoogleTest)
#gtest_discover_tests(myproject)

# Simulation code
add_executable(myproject ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(myproject PUBLIC Eigen3::Eigen gtest gtest_main)
add_test(NAME myproject COMMAND myproject)
