# CMake project file for Fortran
cmake_minimum_required (VERSION 2.8)

enable_language(Fortran)

file(GLOB_RECURSE COMMON_SOURCE *.F90)

add_library(common STATIC ${COMMON_SOURCE})
target_link_libraries(common hdf5_wrapper)
target_link_libraries(common mt_stream)
target_link_libraries(common slps)

target_link_libraries(common ${FFTW3_LIBRARY})
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  target_link_libraries(common iomp5)
  target_link_libraries(common pthread)
  target_link_libraries(common m)
  target_link_libraries(common dl)
endif()

include_directories("${CMAKE_BINARY_DIR}/src/hdf5_wrapper")
include_directories("${CMAKE_BINARY_DIR}/src/mt_stream")
include_directories("${CMAKE_BINARY_DIR}/src/slps")

message(STATUS "BINARY = ${CMAKE_BINARY_DIR}")
