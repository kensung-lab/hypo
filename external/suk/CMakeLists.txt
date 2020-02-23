cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(suk LANGUAGES CXX VERSION 1.0.0)

include(GNUInstallDirs)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug mode")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)  # Compiler-specific features are not used.

# Output directories for binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Options
option(suk_build_executable "Build suk standalone tool" ON)
option(optimise_for_native "Build with -march=native" OFF)

if(optimise_for_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-field-initializers -Wno-missing-braces -pedantic -O3 -march=native")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-field-initializers -Wno-missing-braces -pedantic -O3")
endif()

# build SUK as a static library by default
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")

# Source files
set(LIBSOURCES 
    ${PROJECT_SOURCE_DIR}/include/suk/SolidKmers.hpp
    ${PROJECT_SOURCE_DIR}/src/SolidKmers.cpp
)

# Check Dependencies
find_package (ZLIB REQUIRED)

if (NOT TARGET sdsl)
    add_subdirectory(external/sdsl-lite EXCLUDE_FROM_ALL)    
endif()
if (NOT TARGET slog)
    add_subdirectory(external/slog EXCLUDE_FROM_ALL)
endif() 
if (NOT TARGET kmc_api)
    add_subdirectory(external/kmc_api EXCLUDE_FROM_ALL)
endif()

# Build Library
add_library(suk ${LIBSOURCES})
target_include_directories(suk PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
set_target_properties(suk
    PROPERTIES
    VERSION ${suk_VERSION}
    SOVERSION ${suk_VERSION})
target_link_libraries(suk sdsl slog kmc_api ${ZLIB_LIBRARIES})

# Create target 'install' for installing libraries.
install(TARGETS suk DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/suk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# configure and install pkg-config file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/suk.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/suk-1.pc
    @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/suk-1.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)



if (suk_build_executable)
    # Build binaries.
    set(BINSOURCES 
        ${PROJECT_SOURCE_DIR}/src/main.cpp
    )
    add_executable(suk_bin ${BINSOURCES})

    target_link_libraries(suk_bin suk)
    set_target_properties(suk_bin PROPERTIES OUTPUT_NAME suk)

    install(TARGETS suk_bin DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()



