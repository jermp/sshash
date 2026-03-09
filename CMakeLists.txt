cmake_minimum_required(VERSION 3.5...4.0)
project(SSHASH)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -mavx2")
  if (SSHASH_USE_ARCH_NATIVE AND NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

  if (SSHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

if (CONDA_BUILD)
  set(CONDA_BUILD TRUE)
else()
  set(CONDA_BUILD FALSE)
endif()

option(SSHASH_BUILD_EXECUTABLES "Build sshash executables" ON)
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Conda build: ${CONDA_BUILD}")
MESSAGE(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS "Compiling for processor: ${CMAKE_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "Compiling with flags:${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Build executables: ${SSHASH_BUILD_EXECUTABLES}")

set(Z_LIB_SOURCES
  external/gz/zip_stream.cpp
)

set(CITYHASH_SOURCES
  external/cityhash/cityhash.cpp
)

set(SSHASH_SOURCES
  src/build.cpp
  src/dictionary.cpp
  src/query.cpp
  src/info.cpp
)

set(SSHASH_INCLUDE_DIRS
  external/pthash/include
  external/pthash/external/bits/include
  external/pthash/external/fastmod
  external/pthash/external/bits/external/essentials/include
  external/pthash/external/xxHash
  external/pthash/external/mm_file/include
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Create a static lib
add_library(sshash_static STATIC
  ${Z_LIB_SOURCES}
  ${CITYHASH_SOURCES}
  ${SSHASH_SOURCES}
)

target_include_directories(sshash_static PUBLIC ${SSHASH_INCLUDE_DIRS})

if(SSHASH_BUILD_EXECUTABLES)
  add_executable(sshash tools/sshash.cpp)
  target_include_directories(sshash PUBLIC ${SSHASH_INCLUDE_DIRS})
  target_link_libraries(sshash
    sshash_static
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

  if (CONDA_BUILD)
    install(TARGETS sshash
      RUNTIME DESTINATION bin
    )
  endif()
endif()
