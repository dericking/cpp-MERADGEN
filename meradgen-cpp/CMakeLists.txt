cmake_minimum_required(VERSION 3.14)
project(meradgen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Include directory for C++ headers
set(MERADGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# MERADGEN C++ library
add_library(meradgen_cpp STATIC
  src/cpp/globals.cpp
  src/cpp/meradgen_core.cpp
  src/cpp/meradgen_simpson.cpp
  src/cpp/fsir.cpp
  src/cpp/meradgen_main.cpp
)
target_include_directories(meradgen_cpp PUBLIC ${MERADGEN_INCLUDE_DIR}
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)
option(MERADGEN_PARITY_TRACE "Dump V/T1/Z grid hex traces" OFF)
if(MERADGEN_PARITY_TRACE)
  target_compile_definitions(meradgen_cpp PRIVATE MERADGEN_PARITY_TRACE)
  target_sources(meradgen_cpp PRIVATE src/cpp/meradgen_parity_trace.cpp)
endif()

# Run driver (mirrors run.f)
add_executable(meradgen_run src/driver/run_main.cpp)
target_link_libraries(meradgen_run PRIVATE meradgen_cpp)
target_include_directories(meradgen_run PRIVATE ${MERADGEN_INCLUDE_DIR})

