################################################################################
# Copyright 2022 Matthew Deutsch
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################
cmake_minimum_required(VERSION 3.0) 
project(mineflow)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

################################################################################

if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /W3 /O2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1z")
else()
    message(FATAL_ERROR "Not yet supported")
endif()

################################################################################

add_executable(mineflow
    mineflow.cpp
)
target_compile_definitions(mineflow 
    PUBLIC 
    MVD_MINEFLOW_EXE
)

add_executable(mineflow_test
    mineflow.cpp
)
target_compile_definitions(mineflow_test
    PUBLIC 
    MVD_MINEFLOW_TESTS
)

if(UNIX)
    add_executable(mineflow_test_gmp
        mineflow.cpp
    )
    target_compile_definitions(mineflow_test_gmp
        PUBLIC
        MVD_USE_GMP
        MVD_MINEFLOW_TESTS
    )
    target_link_libraries(mineflow_test_gmp
        gmp
    )
endif()

