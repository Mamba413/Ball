cmake_minimum_required(VERSION 3.12)

include_directories(../lib/gtest-1.7.0/include ../lib/gtest-1.7.0 ../../../src)

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG -fopenmp -O2 -Wall -std=gnu99 -mtune=generic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

add_executable(runBallTest statTest.cpp ../../BD.c ../../BI.c ../../kbcov.c ../../BD.h
        ../../utilities.h ../../utilities.c ../../utilize_cross.h ../../utilize_cross.c ../../Ball_omp.h
        ../../bcor.h ../../bcor.c ../../kbcov.c ../../kbcov.h)

target_link_libraries(runBallTest gtest gtest_main)
