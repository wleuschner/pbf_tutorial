cmake_minimum_required(VERSION 3.5)

project(PSF2D LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(PSF2D
  main.cpp
  mainwindow.cpp
  mainwindow.h
  mainwindow.ui
  canvas.cpp
  canvas.h
  particle.h
  abstractsolver.cpp
  abstractsolver.h
  solver.cpp
  solver.h
  solver_gpu.cpp
  solver_gpu.h
  boundary.h
  abstractkernel.cpp
  abstractkernel.h
  poly6kernel.cpp
  poly6kernel.h
  spikykernel.cpp
  spikykernel.h
  visckernel.cpp
  visckernel.h
  abstractspatialstruct.cpp
  abstractspatialstruct.h
  linearsearch.cpp
  linearsearch.h
  radixsort.cpp
  radixsort.h
  simulation.h
)

find_package(OpenCL)
find_package(OpenMP)

include_directories(${OpenCL_INCLUDE_DIRS})
link_directories(${OpenCL_LIBRARY})
target_link_libraries(PSF2D PRIVATE Qt5::Widgets OpenMP::OpenMP_CXX ${OpenCL_LIBRARY})
