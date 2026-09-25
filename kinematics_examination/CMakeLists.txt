cmake_minimum_required(VERSION 3.14)
project(kinematics_examination LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Prefer meradgen-cpp-final-dev (Approach A / weight_at worktree) when present.
set(_merad_final2 "${CMAKE_CURRENT_SOURCE_DIR}/../meradgen-cpp-final-dev")
set(_merad_final  "${CMAKE_CURRENT_SOURCE_DIR}/../meradgen-cpp-final")
if(EXISTS "${_merad_final2}/CMakeLists.txt")
  set(MERADGEN_TREE "${_merad_final2}")
  message(STATUS "kinematics_examination: linking meradgen-cpp-final-dev")
elseif(EXISTS "${_merad_final}/CMakeLists.txt")
  set(MERADGEN_TREE "${_merad_final}")
  message(STATUS "kinematics_examination: linking meradgen-cpp-final (final-dev missing)")
else()
  message(FATAL_ERROR "No meradgen-cpp-final[-2] tree found next to kinematics_examination/")
endif()

add_subdirectory(${MERADGEN_TREE} ${CMAKE_CURRENT_BINARY_DIR}/meradgen_cpp)

# Detect weight_at without editing the library (another agent owns that API).
set(MERADGEN_HAS_WEIGHT_AT 0)
foreach(_hdr meradgen_molpol.hpp meradgen_api.hpp)
  set(_path "${MERADGEN_TREE}/include/${_hdr}")
  if(EXISTS "${_path}")
    file(READ "${_path}" _hdr_txt)
    string(FIND "${_hdr_txt}" "weight_at" _wa_pos)
    if(NOT _wa_pos EQUAL -1)
      set(MERADGEN_HAS_WEIGHT_AT 1)
    endif()
  endif()
endforeach()
if(MERADGEN_HAS_WEIGHT_AT)
  message(STATUS "kinematics_examination: weight_at DETECTED — full closure enabled")
else()
  message(STATUS "kinematics_examination: weight_at NOT present — scaffold + paired bias only")
endif()

function(ke_target name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE meradgen_cpp)
  target_compile_options(${name} PRIVATE -O2 -ffp-contract=off -fno-fast-math)
  target_compile_definitions(${name} PRIVATE
    MERADGEN_HAS_WEIGHT_AT=${MERADGEN_HAS_WEIGHT_AT})
endfunction()

ke_target(compare_helicity compare_helicity.cpp)
ke_target(closure_approach_a closure_approach_a.cpp)
