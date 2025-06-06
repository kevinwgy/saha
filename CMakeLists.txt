cmake_minimum_required(VERSION 3.10)

# set the project name
project(SahaSolver VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

find_package(Boost 1.72 REQUIRED)
if(NOT "${Boost_INCLUDE_DIR}" STREQUAL "")
  include_directories("${Boost_INCLUDE_DIR}")
endif()

# compiler flags (turn on most of the warnings, except unused variables and parameters
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
include_directories(${CMAKE_SOURCE_DIR}/MathTools)
link_directories(${CMAKE_SOURCE_DIR}/parser)

# compile the parser
add_custom_target(
   extern_lib
   COMMAND $(MAKE)
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/parser
)

# define the macro that allows 1D Riemann solver to output solution
add_compile_definitions(DEBUG_SAHA_SOLVER=1)

# add the executable
add_executable(saha
Main.cpp
IoData.cpp
SahaEquationSolver.cpp
NonIdealSahaEquationSolver.cpp
AtomicIonizationData.cpp
Utils.cpp
MathTools/polynomial_equations.cpp)

# link to libraries
target_link_libraries(saha parser)
##target_link_libraries(m2c petsc mpi parser)
add_dependencies(saha extern_lib)
