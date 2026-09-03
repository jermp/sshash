cmake_minimum_required(VERSION 3.5...4.0)
project(SSHASH)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Options
option(SSHASH_BUILD_EXECUTABLES "Build sshash executables" ON)
option(SSHASH_USE_ARCH_NATIVE "Build with -march=native" OFF)
option(SSHASH_USE_MAX_KMER_LENGTH_63 "Use max kmer length 63" OFF)
option(SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING "Use traditional nucleotide encoding" OFF)
option(SSHASH_USE_SANITIZERS "Build with sanitizers" OFF)

if (UNIX AND (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -mavx2")
    if (SSHASH_USE_ARCH_NATIVE AND NOT CMAKE_CROSSCOMPILING)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif ()
endif ()

if (SSHASH_USE_MAX_KMER_LENGTH_63)
    message(STATUS "SSHash uses a maximum kmer length of 63")
    add_compile_definitions(SSHASH_USE_MAX_KMER_LENGTH_63)
else ()
    message(STATUS "SSHash uses a maximum kmer length of 31")
endif ()

if (SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING)
    message(STATUS "SSHash maps {'A','a'}->0, {'C','c'}->1, {'G','g'}->2, and {'T','t'}->3")
    add_compile_definitions(SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING)
else ()
    message(STATUS "SSHash maps {'A','a'}->0, {'C','c'}->1, {'T','t'}->2, and {'G','g'}->3")
endif ()

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wall -Wextra -Werror -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function -pthread")

    if (SSHASH_USE_SANITIZERS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    endif ()
endif ()

if (CONDA_BUILD)
    set(CONDA_BUILD TRUE)
else ()
    set(CONDA_BUILD FALSE)
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Conda build: ${CONDA_BUILD}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Compiling for processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Compiling with flags:${CMAKE_CXX_FLAGS}")
message(STATUS "Build executables: ${SSHASH_BUILD_EXECUTABLES}")

# Find dependencies
find_package(ZLIB REQUIRED)

# SSHash static library target
add_library(sshash_static STATIC
        external/gz/zip_stream.cpp
        external/cityhash/cityhash.cpp
        src/builder/build.cpp
        src/dictionary.cpp
        src/query.cpp
        src/info.cpp
)

# PUBLIC: Only expose SSHash's own API headers to consuming projects
target_include_directories(sshash_static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/external/bits/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/external/fastmod>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/external/bits/external/essentials/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/external/xxHash>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/pthash/external/mm_file/include>
)

# Link dependencies
target_link_libraries(sshash_static PUBLIC ZLIB::ZLIB)

# Create standard target alias
add_library(sshash::sshash ALIAS sshash_static)

if (SSHASH_BUILD_EXECUTABLES)
    add_executable(sshash tools/sshash.cpp)
    target_link_libraries(sshash PRIVATE sshash_static)

    add_executable(test_alphabet test/test_alphabet.cpp)
    target_link_libraries(test_alphabet PRIVATE sshash_static)

    add_executable(test_minimizer test/test_minimizer.cpp)
    target_link_libraries(test_minimizer PRIVATE sshash_static)

    add_executable(check test/check.cpp)
    target_link_libraries(check PRIVATE sshash_static)

    if (CONDA_BUILD)
        install(TARGETS sshash RUNTIME DESTINATION bin)
    endif ()
endif ()
