cmake_minimum_required(VERSION 3.12)

project(COVID19 LANGUAGES CXX)
set(EXTERN_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(SRC_DIR       "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CODE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/src/epidemics")
set(CUSTOM_CXXFLAGS -Wall -Wextra -march=native)

add_subdirectory(${EXTERN_DIR}/pybind11 sir_pybind11)

pybind11_add_module(
    libepidemics
    ${CODE_DIR}/models/cantons/data.cpp
    ${CODE_DIR}/models/cantons/seiin.cpp
    ${CODE_DIR}/models/cantons/seii_c.cpp
    ${CODE_DIR}/models/cantons/sei_c.cpp
    ${CODE_DIR}/models/cantons/seiin_interventions.cpp
    ${CODE_DIR}/utils/assert.cpp
    ${CODE_DIR}/utils/signal.cpp
    ${CODE_DIR}/bindings/bindings.cpp
    ${CODE_DIR}/bindings/cantons.cpp
)
target_compile_options(libepidemics PUBLIC ${CUSTOM_CXXFLAGS})
target_include_directories(libepidemics PUBLIC ${SRC_DIR})

# Make Release mode the default.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
