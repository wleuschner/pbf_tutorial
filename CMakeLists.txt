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
  solver.cpp
  solver.h
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
)

find_package(OpenMP)
target_link_libraries(PSF2D PRIVATE Qt5::Widgets OpenMP::OpenMP_CXX)
