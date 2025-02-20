
cmake_minimum_required(VERSION 3.0)

project(BrownianMotion)

# ------------------------------------------------------------------------------
# 1) 设置C++标准和编译类型
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 这里可以显式指定Debug或Release，如需可改:
set(CMAKE_BUILD_TYPE Release)

# ------------------------------------------------------------------------------
# 2) 如果有clang可用，则设置为默认编译器
# ------------------------------------------------------------------------------
if(EXISTS /usr/bin/clang AND EXISTS /usr/bin/clang++)
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
endif()

# ------------------------------------------------------------------------------
# 3) 找到pybind11和Python开发头文件
# ------------------------------------------------------------------------------
set(PYBIND11_FINDPYTHON True)
find_package(pybind11 REQUIRED)
include_directories(${pybind11_INCLUDE_DIRS})

find_package(Python3 3.13 REQUIRED COMPONENTS Development)
include_directories(${Python3_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# 4) 定义Python扩展模块 (通过pybind11_add_module)
# ------------------------------------------------------------------------------
# 这里将生成一个名为 "BrownianMotion" 的Python模块
# 最终文件名可能是 BrownianMotion.cpython-3x-darwin.so 或 .pyd (Windows)
pybind11_add_module(BrownianMotion ODModule.cpp)

# 针对Python模块的编译选项
target_compile_options(BrownianMotion PRIVATE -Wall -fPIC)

# 根据CMAKE_BUILD_TYPE区分编译器和链接选项
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(BrownianMotion PRIVATE -fsanitize=address)
    target_link_options(BrownianMotion PRIVATE -fsanitize=address -undefined dynamic_lookup)
else()
    target_compile_options(BrownianMotion PRIVATE -O3)
    target_link_options(BrownianMotion PRIVATE -O3 -undefined dynamic_lookup)
endif()

# ------------------------------------------------------------------------------
# 5) 定义可执行文件 (编译 main.cpp)
# ------------------------------------------------------------------------------
add_executable(BrownianMotionMain main.cpp)

# 针对可执行文件的编译选项
target_compile_options(BrownianMotionMain PRIVATE -Wall -fPIC)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(BrownianMotionMain PRIVATE -fsanitize=address)
    target_link_options(BrownianMotionMain PRIVATE -fsanitize=address)
else()
    target_compile_options(BrownianMotionMain PRIVATE -O3)
    target_link_options(BrownianMotionMain PRIVATE -O3)
endif()

# cmake_minimum_required(VERSION 3.0)

# project(BrownianMotion)

# # set C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CPLUS_INCLUDE_PATH /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1)
# # set default build type to Debug
# set(CMAKE_BUILD_TYPE Release)

# # set compiler to clang, clang++ if available
# if(EXISTS /usr/bin/clang AND EXISTS /usr/bin/clang++)
#     set(CMAKE_C_COMPILER /usr/bin/clang)
#     set(CMAKE_CXX_COMPILER /usr/bin/clang++)
# endif()

# set(PYBIND11_FINDPYTHON True)

# find_package(pybind11 REQUIRED)
# include_directories(${pybind11_INCLUDE_DIRS})

# find_package(Python3 3.13 REQUIRED COMPONENTS Development)
# include_directories(${Python3_INCLUDE_DIRS})

# pybind11_add_module(BrownianMotion ODModule.cpp)

# # add -Wall option
# target_compile_options(BrownianMotion PRIVATE -Wall -fPIC)

# # set fsanitize=address for Debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     target_compile_options(BrownianMotion PRIVATE -fsanitize=address -libstdc++)
#     target_link_options(BrownianMotion PRIVATE -fsanitize=address -undefined dynamic_lookup)
# else()
#     target_compile_options(BrownianMotion PRIVATE -O3)
#     target_link_options(BrownianMotion PRIVATE -O3 -undefined dynamic_lookup)
# endif()

# # add executable
# # add_executable(BrownianMotion main.cpp)