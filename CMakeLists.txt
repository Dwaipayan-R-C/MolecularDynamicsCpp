cmake_minimum_required(VERSION 3.11)

project(myproject)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
set(MYPROJECT_HDRS
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
        ${CMAKE_CURRENT_SOURCE_DIR}/header/domain.h
)

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
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

add_subdirectory(tests)

# Simulation code
add_executable(myproject ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(myproject PUBLIC Eigen3::Eigen MPI::MPI_CXX)
