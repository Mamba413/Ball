cmake_minimum_required(VERSION 3.12)
project(src C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG -fopenmp -O2 -Wall -std=gnu99 -mtune=generic")

add_executable(src test2.c BD.c BI.c utilities.c bcor.c kbcov.c utilize_R.c utilize_cross.c Ball_omp.h)

