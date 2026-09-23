cmake_minimum_required(VERSION 3.14)
project(meradgen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library for Geant4/MolPol. Example binaries are not the deliverable.
option(MERADGEN_BUILD_TOOLS "Build meradgen_run / meradgen_quads (not needed in Geant4)" OFF)

# Only force a build type when this tree is the top-level project.
# As add_subdirectory() from MolPol, inherit Geant4's flags.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

add_library(meradgen_cpp STATIC
  src/cpp/globals.cpp
  src/cpp/meradgen_core.cpp
  src/cpp/meradgen_simpson.cpp
  src/cpp/fsir.cpp
  src/cpp/meradgen_main.cpp
  src/cpp/meradgen_molpol.cpp
)
target_include_directories(meradgen_cpp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
# Keep contraction off even when Geant4 uses -O3 / default fp-contract.
target_compile_options(meradgen_cpp PRIVATE
  -ffp-contract=off -fno-fast-math)

if(MERADGEN_BUILD_TOOLS)
  add_executable(meradgen_run src/driver/run_main.cpp)
  target_link_libraries(meradgen_run PRIVATE meradgen_cpp)
  target_compile_options(meradgen_run PRIVATE
    -O2 -ffp-contract=off -fno-fast-math)

  add_executable(meradgen_quads src/driver/driver_quads.cpp)
  target_link_libraries(meradgen_quads PRIVATE meradgen_cpp)
  target_compile_options(meradgen_quads PRIVATE
    -O2 -ffp-contract=off -fno-fast-math)
endif()
