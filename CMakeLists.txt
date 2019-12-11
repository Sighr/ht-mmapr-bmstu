cmake_minimum_required(VERSION 3.9)
project(mmapr_second_try)

set(CMAKE_CXX_STANDARD 14)

add_executable(mmapr_second_try main.cpp LAE_Solver.cpp Method.cpp)

add_executable(testing Test.cpp LAE_Solver.cpp)