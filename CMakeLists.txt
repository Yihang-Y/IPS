
cmake_minimum_required(VERSION 3.3)

project(BrownianMotion)

# ------------------------------------------------------------------------------
# 1) C++ standard, build type, and compiler
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

# set the include directory globally
include_directories(${CMAKE_SOURCE_DIR}/include)

# ------------------------------------------------------------------------------
# 2) use clang and clang++ if available
# ------------------------------------------------------------------------------
if(EXISTS /usr/bin/clang AND EXISTS /usr/bin/clang++) # what happens if these don't exist??
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
endif()

message("Using C++ compiler: ${CMAKE_CXX_COMPILER}\n")

#  use pybind11?
if(USE_PYBIND)
    message("USE_PYBIND11 is set")

    # ------------------------------------------------------------------------------
    # 3) find pybind11 and Python3
    # ------------------------------------------------------------------------------
    set(PYBIND11_FINDPYTHON True)
    find_package(pybind11 REQUIRED)
    include_directories(${pybind11_INCLUDE_DIRS})
    message("Found pybind11: ${pybind11_INCLUDE_DIRS}")

    find_package(Python3 3.8 REQUIRED COMPONENTS Development)
    include_directories(${Python3_INCLUDE_DIRS})

    # ------------------------------------------------------------------------------
    # 4) add Python module (compile ODModule.cpp)
    # ------------------------------------------------------------------------------
    file(GLOB PYMODULES "pybind/*.cpp")
    foreach(PYMODULE ${PYMODULES})
        get_filename_component(PYMODULE_NAME ${PYMODULE} NAME_WE)
        pybind11_add_module(${PYMODULE_NAME} ${PYMODULE})
        message("Adding Python module: ${PYMODULE_NAME}")
        target_compile_definitions(${PYMODULE_NAME} PRIVATE USE_PYBIND11)
        target_compile_options(${PYMODULE_NAME} PRIVATE -Wall -fPIC)

        # set fsanitize=address for Debug mode
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_options(${PYMODULE_NAME} PRIVATE -fsanitize=address)
            target_link_options(${PYMODULE_NAME} PRIVATE -fsanitize=address -undefined dynamic_lookup)
        else()
            target_compile_options(${PYMODULE_NAME} PRIVATE -O3 -march=native -mtune=native)
            target_link_options(${PYMODULE_NAME} PRIVATE -O3 -undefined dynamic_lookup)
        endif()
        
    endforeach()
endif()

# ------------------------------------------------------------------------------
# 5) add executable (compile all cpp files in exmaple folder)
# ------------------------------------------------------------------------------

file(GLOB EXAMPLES "examples/*.cpp")
foreach(EXAMPLE ${EXAMPLES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE})
    target_compile_options(${EXAMPLE_NAME} PRIVATE -Wall -pthread)

    # set fsanitize=address for Debug mode
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${EXAMPLE_NAME} PRIVATE -fsanitize=address -pthread)
        target_link_options(${EXAMPLE_NAME} PRIVATE -fsanitize=address -pthread)
    else()
        target_compile_options(${EXAMPLE_NAME} PRIVATE -O3 -pthread -march=native -mtune=native)
        target_link_options(${EXAMPLE_NAME} PRIVATE -O3 -pthread)
    endif()
endforeach()
