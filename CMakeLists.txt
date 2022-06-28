cmake_minimum_required(VERSION 3.5)
project(SSHASH)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")        # for multi-threading

  if (SSHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

set(Z_LIB_SOURCES
  include/gz/zip_stream.cpp
)

set(SSHASH_SOURCES
  include/dictionary.cpp
  include/info.cpp
)

# Create a static lib
add_library(sshash_static STATIC
  ${Z_LIB_SOURCES}
  ${SSHASH_SOURCES}
  include/builder/build.cpp
)

# build with the lib
add_executable(build src/build.cpp)
target_link_libraries(build
  sshash_static
  z
)

add_executable(bench src/bench.cpp ${SSHASH_SOURCES})

add_executable(check src/check.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
target_link_libraries(check z)

add_executable(query src/query.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
target_link_libraries(query z)

add_executable(permute src/permute.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
target_link_libraries(permute z)