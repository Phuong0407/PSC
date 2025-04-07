# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(QuadrilateralProject LANGUAGES Fortran)

# Add the quadrilateral module as an object file
add_library(quadrilateral OBJECT quadrilateral.f90)

# Add the executable by linking the test_quadrilateral program with quadrilateral.o
add_executable(test_quadrilateral test_quadrilateral.f90)
target_link_libraries(test_quadrilateral PRIVATE quadrilateral)
