all:
	g++ src_TDEM/main.cc -Wall  -std=c++0x -O3 src_TDEM/allocator.cc src_TDEM/utils.cc src_TDEM/itemGraph.cc src_TDEM/anyoption.cc src_TDEM/sfmt/SFMT.c  -o main_TDEM
