cmake_minimum_required(VERSION 3.6)
project(TMC_OCTv3)

find_package(GSL REQUIRED)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    build/tmcoct_go.o
    build/tmcoct_io.o
    build/tmcoct_mesh.o
    build/tmcoct_nr.o
    build/tmcoct_bias.o
    input/Ellipsoid+2Spheres/BiasStartup.mci
    input/Ellipsoid+2Spheres/Ellipsoid+2Spheres-60Degree.mesh
    input/Ellipsoid+2Spheres/Ellipsoid+2Spheres.opt
    src/tmcoct.h
    src/tmcoct_go.c
    src/tmcoct_io.c
    src/tmcoct_main.c
    src/tmcoct_mesh.c
    src/tmcoct_nr.c
    src/tmcoct_bias.c
    utils/python_utils/mesh_generator.py
    utils/python_utils/requirements.txt
    utils/python_utils/visualization_AScan.py
    utils/python_utils/visualization_BScan.py
    Makefile)

add_executable(TMC_OCTv3 ${SOURCE_FILES})
target_link_libraries(TMC_OCTv3 m)
target_link_libraries(TMC_OCTv3 gsl)
target_link_libraries(TMC_OCTv3 gslcblas)