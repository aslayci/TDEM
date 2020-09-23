cmake_minimum_required(VERSION 3.15)
project(TDEM)

set(CMAKE_CXX_STANDARD 14)

include_directories(src_TDEM)
include_directories(src_TDEM/sfmt)
include_directories(example_input)

add_executable(TDEM
        src_TDEM/sfmt/SFMT-alti.h
        src_TDEM/sfmt/SFMT-common.h
        src_TDEM/sfmt/SFMT-params.h
        src_TDEM/sfmt/SFMT-params11213.h
        src_TDEM/sfmt/SFMT-params1279.h
        src_TDEM/sfmt/SFMT-params132049.h
        src_TDEM/sfmt/SFMT-params19937.h
        src_TDEM/sfmt/SFMT-params216091.h
        src_TDEM/sfmt/SFMT-params2281.h
        src_TDEM/sfmt/SFMT-params4253.h
        src_TDEM/sfmt/SFMT-params44497.h
        src_TDEM/sfmt/SFMT-params607.h
        src_TDEM/sfmt/SFMT-params86243.h
        src_TDEM/sfmt/SFMT-sse2.h
        src_TDEM/sfmt/SFMT.c
        src_TDEM/sfmt/SFMT.h
        src_TDEM/allocator.cc
        src_TDEM/allocator.h
        src_TDEM/anyoption.cc
        src_TDEM/anyoption.h
        src_TDEM/head.h
        src_TDEM/iheap.h
        src_TDEM/itemGraph.cc
        src_TDEM/itemGraph.h
        src_TDEM/main.cc
        src_TDEM/utils.cc
        src_TDEM/utils.h
        example_input/graph.txt
        )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/project/objective-allocation-of-political-ads/codes/TDEM")