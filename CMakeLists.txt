
cmake_minimum_required(VERSION 3.3)

project(BrownianMotion)

# ------------------------------------------------------------------------------
# 1) C++ standard, build type, and compiler
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

# ------------------------------------------------------------------------------
# 2) use clang and clang++ if available
# ------------------------------------------------------------------------------
if(EXISTS /usr/bin/clang AND EXISTS /usr/bin/clang++)
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
endif()

# ------------------------------------------------------------------------------
# 3) find pybind11 and Python3
# ------------------------------------------------------------------------------
set(PYBIND11_FINDPYTHON True)
find_package(pybind11 REQUIRED)
include_directories(${pybind11_INCLUDE_DIRS})

find_package(Python3 3.13 REQUIRED COMPONENTS Development)
include_directories(${Python3_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# 4) add Python module (compile ODModule.cpp)
# ------------------------------------------------------------------------------
pybind11_add_module(BrownianMotion ODModule.cpp)

# compile options for the Python module
target_compile_options(BrownianMotion PRIVATE -Wall -fPIC)

# set fsanitize=address for Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(BrownianMotion PRIVATE -fsanitize=address)
    target_link_options(BrownianMotion PRIVATE -fsanitize=address -undefined dynamic_lookup)
else()
    target_compile_options(BrownianMotion PRIVATE -O3)
    target_link_options(BrownianMotion PRIVATE -O3 -undefined dynamic_lookup)
endif()

# ------------------------------------------------------------------------------
# 5) add executable (compile main.cpp)
# ------------------------------------------------------------------------------
add_executable(BrownianMotionMain main.cpp)

# compile options for the executable
target_compile_options(BrownianMotionMain PRIVATE -Wall -fPIC)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(BrownianMotionMain PRIVATE -fsanitize=address)
    target_link_options(BrownianMotionMain PRIVATE -fsanitize=address -undefined dynamic_lookup)
else()
    target_compile_options(BrownianMotionMain PRIVATE -O3)
    target_link_options(BrownianMotionMain PRIVATE -O3 -undefined dynamic_lookup)
endif()