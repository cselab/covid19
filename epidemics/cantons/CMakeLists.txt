cmake_minimum_required(VERSION 3.12)

# Get compiler, cxxflags and linker flags from Korali.
execute_process(
    COMMAND python3 -m korali.cxx --compiler
    OUTPUT_VARIABLE KORALI_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
    COMMAND python3 -m korali.cxx --cflags
    OUTPUT_VARIABLE KORALI_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
    COMMAND python3 -m korali.cxx --libs
    OUTPUT_VARIABLE KORALI_LINKERFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
# Remove the erronous -Wl,.../build/lib from the linker flags.
string(REGEX REPLACE " -Wl,[^ ]*build/lib" "" KORALI_LINKERFLAGS ${KORALI_LINKERFLAGS})
# String to list.
string(REPLACE " " ";" KORALI_CXXFLAGS ${KORALI_CXXFLAGS})
string(REPLACE " " ";" KORALI_LINKERFLAGS ${KORALI_LINKERFLAGS})

# Compiler seems to be able configurable globally.
# This has to be done BEFORE project().
set(CMAKE_CXX_COMPILER "${KORALI_COMPILER}")

project(COVID19 LANGUAGES CXX)
set(EXTERN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CUSTOM_CXXFLAGS -Wall -Wextra -march=native)

# Files common to .so and to the executable.
SET(COMMON_SRC_FILES
    ${SRC_DIR}/common.cpp
    ${SRC_DIR}/data.cpp
    ${SRC_DIR}/model_seiin.cpp
    ${SRC_DIR}/model_seii_c.cpp
    ${SRC_DIR}/model_sei_c.cpp
    ${SRC_DIR}/utils.cpp)

# Set up the Python .so code.
add_subdirectory(${EXTERN_DIR}/pybind11 sir_pybind11)
pybind11_add_module(
    libsolver
    ${COMMON_SRC_FILES}
    ${SRC_DIR}/bindings.cpp)
target_compile_options(libsolver PUBLIC ${CUSTOM_CXXFLAGS})

# Set up the C++-only code.
add_executable(
    solver
    ${COMMON_SRC_FILES}
    ${SRC_DIR}/korali_setup.cpp)

target_compile_options(solver PUBLIC ${KORALI_CXXFLAGS})
target_compile_options(solver PUBLIC ${CUSTOM_CXXFLAGS} -Wno-unused-parameter)
set_target_properties(solver PROPERTIES LINK_FLAGS ${KORALI_LINKERFLAGS})

# Make Release mode the default.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
