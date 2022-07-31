cmake_minimum_required(VERSION 3.11)
set( CMAKE_CXX_COMPILER "/usr/bin/g++" )
set( CMAKE_C_COMPILER "/usr/bin/gcc" )
project(yamd)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)
# #Adding gtset
# include(FetchContent)
# FetchContent_Declare(
#         googletest
#         URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
# )
# FetchContent_MakeAvailable(googletest)
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#////////
add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(MYPROJECT_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/header/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/atom_structure.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/xyz.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/force.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/lj_direct_summation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/lj_kinetic_and_potential.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/berendsen_thermostat.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/gupta.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/mpi_support.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/variable_def.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/neighbors.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/milestones.h

        )

#        )

set(MYPROJECT_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/verlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/berendsen_thermostat.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/domain.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/gupta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/lj_direct_summation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/lj_kinetic_potential.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/milestones.cpp
        )

set(MYPROJECT_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)


# Simulation code
add_executable(myproject ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(myproject PUBLIC Eigen3::Eigen MPI::MPI_CXX)
# MPI::MPI_CXX

