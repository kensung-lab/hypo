cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(slog LANGUAGES CXX VERSION 1.0.0)

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -Wextra -pedantic ")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Compiler-specific features are not used.

# Output directories for binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Source files
set(SOURCES 
    "${PROJECT_SOURCE_DIR}/src/Monitor.cpp"
    "${PROJECT_SOURCE_DIR}/src/memory_usage.h"
)

# Build library.
add_library(slog STATIC ${SOURCES})
# Enable .cpp files to add headers in the include folder.
include_directories (
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)
target_include_directories(slog PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    install(TARGETS slog DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/slog DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    


