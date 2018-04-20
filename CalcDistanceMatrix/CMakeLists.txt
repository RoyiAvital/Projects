cmake_minimum_required (VERSION 3.10)

# General Settings
set(CMAKE_VERBOSE_MAKEFILE ON)

# Project
project(CalcDistanceMatrix)

set(EIGEN_INCLUDE_FOLDER "${CMAKE_SOURCE_DIR}/../../../Libraries/EigenLib/Eigen/")

set(HEADER_FILES "$EIGEN_INCLUDE_FOLDER/Dense")

set(RESOURCES "")

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/CalcDistanceMatrix/*.cpp")
list(APPEND SOURCES_EXE ${RESOURCES} ${SOURCES})

# Removing items for DLL / LIB See https://stackoverflow.com/a/16450083
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/CalcDistanceMatrix/CalcDistanceMatrixMain.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/CalcDistanceMatrix/CalcDistanceMatrixUnitTest.cpp")
list(APPEND SOURCES_LIB ${RESOURCES} ${SOURCES})

message("SOURCES LIB: ${SOURCES_LIB}")

# Targets
# add_executable(CalcDistanceMatrix ${SOURCES_EXE} ${HEADER_FILES})
# add_library(CalcDistanceMatrixDll SHARED ${SOURCES_LIB} ${HEADER_FILES})
add_executable(CalcDistanceMatrix ${SOURCES_EXE})
add_library(CalcDistanceMatrixDll SHARED ${SOURCES_LIB})

# Set Compiling Flags
target_compile_options(CalcDistanceMatrix PUBLIC -Ofast -fopenmp -mavx2)
target_compile_definitions(CalcDistanceMatrix PUBLIC -D_DNDEBUG)
target_compile_options(CalcDistanceMatrixDll PUBLIC -Ofast -fopenmp -mavx2)
target_compile_definitions(CalcDistanceMatrixDll PUBLIC -D_USRDLL -D_DNDEBUG)

# For LINK_FLAGS use single "...."
SET_TARGET_PROPERTIES(CalcDistanceMatrix PROPERTIES LINK_FLAGS -fopenmp) # Set Link Flags
SET_TARGET_PROPERTIES(CalcDistanceMatrixDll PROPERTIES LINK_FLAGS -fopenmp) # Set Link Flags

SET_TARGET_PROPERTIES(CalcDistanceMatrixDll PROPERTIES PREFIX "") # Remove lib Prefix
SET_TARGET_PROPERTIES(CalcDistanceMatrix PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/GCC/") # Remove lib Prefix
SET_TARGET_PROPERTIES(CalcDistanceMatrixDll PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/GCC/") # DLL is considered Run Time
# Currently works with GCC only

# Run the Process
# cmake .. -G"Visual Studio 15 2017 Win64"
# cmake .. -G"Visual Studio 15 2017 Win64" -T"Intel C++ Compiler 18.0"
# cmake .. -G"MinGW Makefiles"

# cmake --build . --config Release



