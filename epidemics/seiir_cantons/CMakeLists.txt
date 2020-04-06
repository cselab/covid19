cmake_minimum_required(VERSION 3.10)

project(SIR VERSION LANGUAGES CXX)
set(EXTERN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_subdirectory(${EXTERN_DIR}/pybind11 sir_pybind11)
pybind11_add_module(
    libseiir
    ${SRC_DIR}/model.cpp
    ${SRC_DIR}/bindings.cpp)
