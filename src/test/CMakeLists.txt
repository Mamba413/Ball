cmake_minimum_required(VERSION 3.12)

project(ball_tests)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG -fopenmp -O2 -Wall -std=gnu99 -mtune=generic")

add_subdirectory(lib/gtest-1.7.0)
add_subdirectory(ball_test)