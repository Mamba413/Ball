cmake_minimum_required(VERSION 3.12)
project(ball C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG -fopenmp -O2 -Wall -std=gnu99 -mtune=generic")

add_executable(src BD.c BI.c utilities.c bcor.c kbcov.c Ball_omp.h kbd.c kbd.h)

add_subdirectory(test)
include_directories(test/lib/gtest/include test/lib/gtest)