cmake_minimum_required(VERSION 3.5...4.0)
project(SSHASH)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  if (SSHASH_USE_ARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2")
endif()

if (SSHASH_USE_MAX_KMER_LENGTH_63)
  MESSAGE(STATUS "SSHash uses a maximum kmer length of 63")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSSHASH_USE_MAX_KMER_LENGTH_63")
else()
  MESSAGE(STATUS "SSHash uses a maximum kmer length of 31")
endif()

if (SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING)
  MESSAGE(STATUS "SSHash maps {'A','a'}->0, {'C','c'}->1, {'G','g'}->2, and {'T','t'}->3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING")
else()
  MESSAGE(STATUS "SSHash maps {'A','a'}->0, {'C','c'}->1, {'T','t'}->2, and {'G','g'}->3")
endif()

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

  if (SSHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Compiling for processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "Compiling with flags:${CMAKE_CXX_FLAGS}")

include_directories(.) # all include paths relative to parent directory
include_directories(external/pthash/include)
include_directories(external/pthash/external/bits/include)
include_directories(external/pthash/external/fastmod)
include_directories(external/pthash/external/bits/external/essentials/include)
include_directories(external/pthash/external/xxHash)
include_directories(external/pthash/external/mm_file/include)

set(Z_LIB_SOURCES
  external/gz/zip_stream.cpp
)

set(SSHASH_SOURCES
  src/build.cpp
  src/dictionary.cpp
  src/query.cpp
  src/dump.cpp
  src/info.cpp
  src/statistics.cpp
)

# Create a static lib
add_library(sshash_static STATIC
  ${Z_LIB_SOURCES}
  ${SSHASH_SOURCES}
)

add_executable(sshash tools/sshash.cpp)
target_link_libraries(sshash
  z
)

# tests:

add_executable(test_alphabet test/test_alphabet.cpp)
target_link_libraries(test_alphabet
  sshash_static
)

add_executable(check test/check.cpp)
target_link_libraries(check
  sshash_static
  z
)
