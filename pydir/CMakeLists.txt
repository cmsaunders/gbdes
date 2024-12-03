cmake_minimum_required(VERSION 3.20)
project(gbdes_py)

find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

set(PYTHON_MODULE_EXTENSION ".so")

pybind11_add_module(wcsfit wcsfit.cpp)

target_link_libraries(wcsfit PUBLIC gbdes)

target_include_directories(wcsfit PUBLIC ${PROJECT_SOURCE_DIR}/include ${ASTSHIM_INCLUDES})

# add import test
add_test(NAME PythonImportTest
         COMMAND ${Python3_EXECUTABLE} -c "import wcsfit" && printenv)

# make sure we don't use the newly build shlib for testing
set_tests_properties(PythonImportTest PROPERTIES
                     WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                     ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pydir")

add_custom_command(TARGET wcsfit POST_BUILD
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running import test after build"
)

install(TARGETS wcsfit DESTINATION "${PROJECT_SOURCE_DIR}")
