cmake_minimum_required(VERSION 3.14)
project(meradgen_parity_harness LANGUAGES C CXX Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 14)

set(PARITY_FFLAGS -O2 -ffp-contract=off -fno-fast-math -fno-automatic
                  -ffixed-line-length-none -fdollar-ok -fno-second-underscore)
set(PARITY_CXXFLAGS -O2 -ffp-contract=off -fno-fast-math)

get_filename_component(MERADGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(FORTRAN_REF "${MERADGEN_ROOT}/meradgen-fortran")
set(CPP_ROOT "${MERADGEN_ROOT}/meradgen-cpp")
set(FORTRAN_COPY "${CMAKE_CURRENT_BINARY_DIR}/fortran_src")

option(PARITY_TRACE "Instrument grids + vacpol/L1f (needs a local testing-directory)" OFF)
set(MERADGEN_PARITY_TRACE ${PARITY_TRACE} CACHE BOOL "Enable meradgen-cpp traces" FORCE)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

if(PARITY_TRACE)
  message(FATAL_ERROR
    "PARITY_TRACE needs an instrumented FORTRAN tree; validation_checks_new "
    "was removed. Rebuild with -DPARITY_TRACE=OFF (default).")
else()
  add_custom_command(
    OUTPUT
      ${FORTRAN_COPY}/meradgen10.f
      ${FORTRAN_COPY}/fsir.f
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/inject_rand4.py
            --src ${FORTRAN_REF}
            --dst ${FORTRAN_COPY}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/inject_rand4.py
      ${FORTRAN_REF}/meradgen10.f
      ${FORTRAN_REF}/fsir.f
    COMMENT "Copy meradgen-fortran and inject rand4 (build tree only)"
  )
  add_custom_target(fortran_rand4_src DEPENDS
    ${FORTRAN_COPY}/meradgen10.f
    ${FORTRAN_COPY}/fsir.f
  )
  set_source_files_properties(
    ${FORTRAN_COPY}/meradgen10.f
    ${FORTRAN_COPY}/fsir.f
    PROPERTIES GENERATED TRUE
  )
  add_executable(meradgen_parity_fortran
    ${FORTRAN_COPY}/meradgen10.f
    ${FORTRAN_COPY}/fsir.f
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/driver_parity.f
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/parity_trace.f
  )
  add_dependencies(meradgen_parity_fortran fortran_rand4_src)
  target_include_directories(meradgen_parity_fortran PRIVATE ${FORTRAN_COPY}/include)
endif()
target_compile_options(meradgen_parity_fortran PRIVATE ${PARITY_FFLAGS})

add_subdirectory(${CPP_ROOT} meradgen_cpp_subbuild)
add_executable(meradgen_parity_cpp cpp/driver_parity.cpp)
target_link_libraries(meradgen_parity_cpp PRIVATE meradgen_cpp)
target_include_directories(meradgen_parity_cpp PRIVATE ${CPP_ROOT}/include)
target_compile_options(meradgen_parity_cpp PRIVATE ${PARITY_CXXFLAGS})
target_compile_options(meradgen_cpp PRIVATE ${PARITY_CXXFLAGS})
