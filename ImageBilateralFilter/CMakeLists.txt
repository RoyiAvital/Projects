cmake_minimum_required (VERSION 3.10)

# General Settings Prior to Project
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

message("Release Default Flags: ${CMAKE_C_FLAGS_RELEASE}")
message("Release Default Definitions: ${COMPILE_DEFINITIONS_RELEASE}")

# Project
project(ImageBilateralFilter)

# General Settings Post to Project
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Defining Paths, Sources and Resources
set(CODE_FILES_PATH "${CMAKE_SOURCE_DIR}/ImageBilateralFilter/")

set(RESOURCES "")

file(GLOB SOURCES "${CODE_FILES_PATH}*.c*")
list(REMOVE_ITEM SOURCES "${CODE_FILES_PATH}GaussianFilterFa.c")


list(APPEND SOURCES_EXE ${RESOURCES} ${SOURCES})

# Removing items for DLL / LIB See https://stackoverflow.com/a/16450083
list(REMOVE_ITEM SOURCES "${CODE_FILES_PATH}ImageBilateralFilterMain.c")
list(REMOVE_ITEM SOURCES "${CODE_FILES_PATH}ImageBilateralFilterUnitTest.c")
list(APPEND SOURCES_LIB ${RESOURCES} ${SOURCES})

# Adding Libraries


# Information
message("\n")

message("SOURCE FILES EXE:")
message("\n")
FOREACH(EXE_FILE ${SOURCES_EXE})
    message("SOURCES EXE: ${EXE_FILE}")
ENDFOREACH(EXE_FILE)
message("\n")

message("SOURCE FILES LIB:")
message("\n")
FOREACH(LIB_FILE ${SOURCES_LIB})
    message("SOURCES LIB: ${LIB_FILE}")
ENDFOREACH(LIB_FILE)
message("\n")

# Set parameters by Compiler
# https://stackoverflow.com/questions/10046114/in-cmake-how-can-i-test-if-the-compiler-is-clang

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    message("Using LLVM CLANG Compiler")
	message("\n")

	SET(OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/CLANG/")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    message("Using GNU GCC Compiler")
	message("\n")
	
	if (WIN32)
		# Windows
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _USRDLL)
	elseif(APPLE)
		# macOS
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _DYLIB)
	elseif(UNIX)
		# Unix / Linux
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _SO)
	endif()
	SET(COMPILER_DEFS_EXE NDEBUG)
	SET(COMPILER_FLAGS -Ofast -fopenmp -mavx2 -static-libgcc)
	# The flag -static-libgcc will also make OpenMP static
	SET(LINKER_FLAGS "-fopenmp -static-libgcc")
	
	SET(OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/GCC/")
	
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    message("Using Intel ICC Compiler")
	message("\n")

	if (WIN32)
		# Windows
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _USRDLL)
		SET(COMPILER_FLAGS -O3 /Qopenmp /arch:avx2)
		SET(LINKER_FLAGS /Qopenmp)
	elseif(APPLE)
		# macOS
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _DYLIB)
		SET(COMPILER_FLAGS -Ofast -qopenmp -mavx2)
		SET(LINKER_FLAGS -qopenmp)

		# ICC has issues with uint16_t on macOS
		# SET_SOURCE_FILES_PROPERTIES(${SOURCES_EXE} PROPERTIES LANGUAGE CXX)
		# SET_SOURCE_FILES_PROPERTIES(${SOURCES_LIB} PROPERTIES LANGUAGE CXX)
	elseif(UNIX)
		# Unix / Linux
		SET(COMPILER_DEFS_SHARED_LIB NDEBUG _SO)
	endif()
	SET(COMPILER_DEFS_EXE NDEBUG)
	
	SET(OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/ICC/")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    message("Using Microsoft MSVC Compiler")
	message("\n")

	SET(OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/Release/")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    # using  AppleClang
    message("Using Apple CLANG Compiler")
	message("\n")

	SET(OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/x64/CLANG/")

endif()

# Link Intel Compiler + Intel MKL
# https://www.scivision.co/intel-compiler-cmake-make-load-libraries/

# Add Batch File Commands
# https://stackoverflow.com/questions/10383021


# Targets
add_executable(ImageBilateralFilterExe ${SOURCES_EXE}) #Executable
add_library(ImageBilateralFilterDll SHARED ${SOURCES_LIB}) # Dynamic Library
add_library(ImageBilateralFilterLib STATIC ${SOURCES_LIB}) # Static Library

# Set Compiling Flags
target_compile_options(ImageBilateralFilterExe PUBLIC ${COMPILER_FLAGS})
target_compile_definitions(ImageBilateralFilterExe PUBLIC ${COMPILER_DEFS_EXE})
target_compile_options(ImageBilateralFilterDll PUBLIC ${COMPILER_FLAGS})
target_compile_definitions(ImageBilateralFilterDll PUBLIC ${COMPILER_DEFS_SHARED_LIB})
target_compile_options(ImageBilateralFilterLib PUBLIC ${COMPILER_FLAGS})
target_compile_definitions(ImageBilateralFilterLib PUBLIC ${COMPILER_DEFS_SHARED_LIB})

# Set Link Flags
SET_TARGET_PROPERTIES(ImageBilateralFilterExe PROPERTIES LINK_FLAGS ${LINKER_FLAGS})
SET_TARGET_PROPERTIES(ImageBilateralFilterDll PROPERTIES LINK_FLAGS ${LINKER_FLAGS})
SET_TARGET_PROPERTIES(ImageBilateralFilterLib PROPERTIES LINK_FLAGS ${LINKER_FLAGS})

# Set Output Name
# See https://stackoverflow.com/questions/31038963
SET_TARGET_PROPERTIES(ImageBilateralFilterExe PROPERTIES OUTPUT_NAME "${CMAKE_PROJECT_NAME}") 
SET_TARGET_PROPERTIES(ImageBilateralFilterDll PROPERTIES PREFIX "") # Remove lib Prefix
SET_TARGET_PROPERTIES(ImageBilateralFilterDll PROPERTIES OUTPUT_NAME "${CMAKE_PROJECT_NAME}")
SET_TARGET_PROPERTIES(ImageBilateralFilterLib PROPERTIES PREFIX "") # Remove lib Prefix
SET_TARGET_PROPERTIES(ImageBilateralFilterLib PROPERTIES OUTPUT_NAME "${CMAKE_PROJECT_NAME}") # See https://stackoverflow.com/questions/31038963

# Set Output Directory
SET_TARGET_PROPERTIES(ImageBilateralFilterExe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY})
SET_TARGET_PROPERTIES(ImageBilateralFilterDll PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}) # DLL is considered Run Time
SET_TARGET_PROPERTIES(ImageBilateralFilterDll PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}) # Library
SET_TARGET_PROPERTIES(ImageBilateralFilterLib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}) # Static Library is considered Archive

# Generation Commands:

# Windows
# MS Build with MSVC Compiler:
# cmake -G"Visual Studio 15 2017 Win64" ..
# MS Build with Intel Compiler:
# cmake -G"Visual Studio 15 2017 Win64" -T"Intel C++ Compiler 18.0" ..
# MinGW with GCC:
# cmake -G"MinGW Makefiles" ..
# Ninja with Intel Compiler (Requires running "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.3.210\windows\bin\ipsxe-comp-vars.bat" intel64 vs2017)
# cmake -G"Ninja" -D CMAKE_C_COMPILER=icl -D CMAKE_CXX_COMPILER=icl ..
# Ninja with CLang (LLVM) Compiler (Requires running the MS Build Environment [x64 Native Tools Command Prompt of VS], See https://stackoverflow.com/questions/46553436):
# cmake -G"Ninja" -D CMAKE_C_COMPILER=clang-cl -D CMAKE_CXX_COMPILER=clang-cl -D CMAKE_LINKER=D:\Applications\Programming\LLVM\bin\clang-cl ..
# cmake -G"Ninja" -D CMAKE_C_COMPILER=D:/Applications/Programming/LLVM/bin/clang-cl.exe -D CMAKE_CXX_COMPILER=D:/Applications/Programming/LLVM/bin/clang-cl.exe -D CMAKE_LINKER=D:/Applications/Programming/LLVM/bin/lld-link.exe ..

# macOS
# Ninja with GCC ():
# cmake -G"Ninja" -D CMAKE_C_COMPILER=gcc-8.2 -D CMAKE_CXX_COMPILER=g++-8.2 ..
# Ninja with Intel Compiler:
# cmake -G"Ninja" -D CMAKE_C_COMPILER=icc -D CMAKE_CXX_COMPILER=icpc ..
# XCode Build with Intel Compiler (Does it work?):
# cmake -G"Xcode" -T"Intel C++ Compiler 18.0"

# Build Command
# cmake --build . --config Release

# To Do
# Replicate Command Line which enables Clang (LLVM) + Ninja

