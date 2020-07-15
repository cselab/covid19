cmake_minimum_required(VERSION 3.12)

project(COVID19 LANGUAGES CXX)

option(ENABLE_BACKWARD_CPP "Whether or not to use backward-cpp, for automatically printing stack traces." OFF)

set(EXTERN_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(SRC_DIR       "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CODE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/src/epidemics")
set(CUSTOM_CXXFLAGS -Wall -Wextra -march=native)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)  # Use -std=c++14, not -std=gnu++14.
add_subdirectory(${EXTERN_DIR}/pybind11 epidemics_pybind11)
if(ENABLE_BACKWARD_CPP)
    add_subdirectory(${EXTERN_DIR}/backward-cpp epidemics_backward_cpp)
endif()

include(GoogleTest)
option(INSTALL_GTEST OFF)
add_subdirectory(${EXTERN_DIR}/googletest)

find_package(Boost)

set(GENERATED_COUNTRY_BINDINGS
    ${CODE_DIR}/bindings/_cantons.generated.cpp
    ${CODE_DIR}/bindings/_cantons.sei_c.generated.cpp
    ${CODE_DIR}/bindings/_cantons.seii_c.generated.cpp
    ${CODE_DIR}/bindings/_cantons.seiin.generated.cpp
    ${CODE_DIR}/bindings/_cantons.seiin_interventions.generated.cpp
    ${CODE_DIR}/bindings/_country.generated.cpp
    ${CODE_DIR}/bindings/_country.sir.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_int.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_int_r0.generated.cpp
    ${CODE_DIR}/bindings/_country.seir.generated.cpp
    ${CODE_DIR}/bindings/_country.seir_int.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir_int.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.sird_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_intexp_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.sir_intrem_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seir_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seir_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seird_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seir_intexp_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seir_intrem_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seiird_int_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir_intexp_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.seiir_intrem_reparam.generated.cpp
    ${CODE_DIR}/bindings/_country.cz_int.generated.cpp

)

add_custom_command(
    OUTPUT
    ${GENERATED_COUNTRY_BINDINGS}
    COMMAND
    ${PYTHON_EXECUTABLE}
    ${CODE_DIR}/bindings/generate_bindings.py
    DEPENDS
    ${CODE_DIR}/bindings/generate_bindings.py
    ${CODE_DIR}/bindings/cantons.template.cpp
    ${CODE_DIR}/bindings/cantons.model.template.cpp
    ${CODE_DIR}/bindings/country.template.cpp
    ${CODE_DIR}/bindings/country.model.template.cpp
    COMMENT "Generating bindings for cantons and country models."
)

add_library(
    libepidemics_core
    STATIC
    # For now, bindings are responsible for instantiating template
    # functions `Solver::solve`, so there's not that many files here.
    ${CODE_DIR}/models/cantons/data.cpp
    ${CODE_DIR}/utils/assert.cpp
    ${CODE_DIR}/utils/backward.cpp
    ${CODE_DIR}/utils/signal.cpp
)
set_property(TARGET libepidemics_core PROPERTY POSITION_INDEPENDENT_CODE ON)
target_compile_options(libepidemics_core PUBLIC ${CUSTOM_CXXFLAGS})
target_include_directories(libepidemics_core PUBLIC ${SRC_DIR})
# Backwards-compatible alias to Boost::headers.
target_link_libraries(libepidemics_core PUBLIC Boost::boost)

add_library(libepidemics::core ALIAS libepidemics_core)

pybind11_add_module(
    libepidemics
    ${GENERATED_COUNTRY_BINDINGS}
    ${CODE_DIR}/bindings/bindings.cpp
)
target_link_libraries(libepidemics PRIVATE libepidemics::core)


###########################
# Test
###########################
set(UNIT_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp")
add_executable(
    libepidemics_unittests
    ${UNIT_SRC_DIR}/test_autodiff.cpp
)

target_link_libraries(libepidemics_unittests PRIVATE gtest gtest_main libepidemics::core)
gtest_discover_tests(libepidemics_unittests)

if(ENABLE_BACKWARD_CPP)
    target_compile_definitions(libepidemics_core PUBLIC EPIDEMICS_ENABLE_BACKWARD_CPP)
    add_backward(libepidemics_core)
    add_backward(libepidemics)
    add_backward(libepidemics_unittests)
endif()

# Make Release mode the default.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
