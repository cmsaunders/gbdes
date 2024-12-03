cmake_minimum_required(VERSION 3.20)
project(gbdes_py)

find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

set(PYTHON_MODULE_EXTENSION ".so")

pybind11_add_module(wcsfit wcsfit.cpp)

target_link_libraries(wcsfit PUBLIC gbdes)

target_include_directories(wcsfit PUBLIC ${PROJECT_SOURCE_DIR}/include ${ASTSHIM_INCLUDES})

install(TARGETS wcsfit DESTINATION "${PROJECT_SOURCE_DIR}")
