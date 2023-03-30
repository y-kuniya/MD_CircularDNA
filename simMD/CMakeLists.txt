cmake_minimum_required(VERSION 3.1)
project(simulationMD CXX)

add_library(simMD STATIC 
        ./src/variables.cpp
        ./src/observer.cpp 
        ./src/tools.cpp
        ./src/force.cpp
        ./src/md.cpp)

target_compile_features(simMD PRIVATE cxx_std_11)
target_include_directories(simMD INTERFACE ./include)