cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(hypo LANGUAGES CXX VERSION 1.0.0)

include(GNUInstallDirs)

#set( CMAKE_VERBOSE_MAKEFILE on )

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
endif()

option(optimise_for_native "Build with -march=native" OFF)

if(optimise_for_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3 -march=native")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3")
endif()



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)  # Compiler-specific features are not used.

# Output directories for binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Source files
set(SOURCES 
    ${PROJECT_SOURCE_DIR}/include/Alignment.hpp
    ${PROJECT_SOURCE_DIR}/include/Contig.hpp
    ${PROJECT_SOURCE_DIR}/include/Filter.hpp
    ${PROJECT_SOURCE_DIR}/include/globalDefs.hpp
    ${PROJECT_SOURCE_DIR}/include/Hypo.hpp
    ${PROJECT_SOURCE_DIR}/include/kseq.h
    ${PROJECT_SOURCE_DIR}/include/MinimizerDeque.hpp
    ${PROJECT_SOURCE_DIR}/include/PackedSeq.hpp
    ${PROJECT_SOURCE_DIR}/include/Window.hpp
    ${PROJECT_SOURCE_DIR}/src/Alignment.cpp
    ${PROJECT_SOURCE_DIR}/src/Contig.cpp
    ${PROJECT_SOURCE_DIR}/src/Contig.cpp
    ${PROJECT_SOURCE_DIR}/src/Hypo.cpp
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/PackedSeq.cpp
    ${PROJECT_SOURCE_DIR}/src/Window.cpp
)

# Build binaries.
add_executable(hypo ${SOURCES})

# Enable .cpp files to add headers in the include folder.
include_directories (
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

# Check Dependencies
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Install External libs

if(CMAKE_BUILD_TYPE MATCHES Conda)
    message("Conda mode")
    find_library(HTS hts)
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/kseq.h HINTS "${HTSLIB_INSTALL_DIR}")
    if (HTS AND HTSLIB_INCLUDE_DIR)
        message(STATUS "   htslib and headers found: ${HTS} and ${HTSLIB_INCLUDE_DIR}")
        include_directories(${HTSLIB_INCLUDE_DIR})
        link_directories(${HTS})
    else()
        message(FATAL_ERROR "   htslib NOT found (in Conda environment). Please try buiding from the source.")
    endif()
    
    find_library(SDSL sdsl)
    find_path(SDSLLIB_INCLUDE_DIR NAMES sdsl/bit_vectors.hpp HINTS "${SDSLLIB_INSTALL_DIR}")
    if (SDSL AND SDSLLIB_INCLUDE_DIR)
        message(STATUS "  sdsllib and headers found: ${SDSL} and ${SDSLLIB_INCLUDE_DIR}")
        include_directories(${HTSLIB_INCLUDE_DIR})
        link_directories(${HTS})
    else()
        message(FATAL_ERROR "   sdsllib NOT found (in Conda environment). Please try buiding from the source.")
    endif()
else()
    message("Installing from Source")
    set(htslib_root ${PROJECT_SOURCE_DIR}/external/install/htslib)
    set(htslib_static_dir ${htslib_root}/lib)
    find_library(HTSFOUND hts PATHS ${htslib_static_dir})
    if (NOT HTSFOUND)
        message(STATUS "   Assuming htslib version 1.10 or higher is installed.")
        find_library(HTS hts)
        find_path(HTSLIB_INCLUDE_DIR NAMES htslib/kseq.h HINTS "${HTSLIB_INSTALL_DIR}")
        if (HTS AND HTSLIB_INCLUDE_DIR)
            message(STATUS "   htslib and headers found: ${HTS} and ${HTSLIB_INCLUDE_DIR}")
            include_directories(${HTSLIB_INCLUDE_DIR})
            link_directories(${HTS})
        else()
            message(FATAL_ERROR "   htslib NOT found (locally or globally). Please run install_deps.sh provided in the hypo folder.")
        endif() 
    else()
        message(STATUS "   Local htslib is being used.")
        add_library(htslib UNKNOWN IMPORTED) # or STATIC instead of SHARED
        set_target_properties(htslib PROPERTIES
            IMPORTED_LOCATION "${htslib_static_dir}/libhts.so"
            INTERFACE_INCLUDE_DIRECTORIES "${htslib_root}/include"
        )   
        set(HTS htslib)
    endif()

    if (NOT TARGET sdsl)
        add_subdirectory(external/sdsl-lite EXCLUDE_FROM_ALL)    
    endif()
endif()

if (NOT TARGET spoa)
    add_subdirectory(external/spoa EXCLUDE_FROM_ALL)
endif()

if (NOT TARGET suk)
    add_subdirectory(external/suk EXCLUDE_FROM_ALL)
endif() 
if (NOT TARGET slog)
    add_subdirectory(external/slog EXCLUDE_FROM_ALL)
endif() 

# This maybe needed for filesystem: target_link_libraries
#target_link_libraries(hypo PRIVATE spoa pthread sdsl suk slog ${ZLIB_LIBRARIES} ${HTS} OpenMP::OpenMP_CXX "stdc++fs")
target_link_libraries(hypo PRIVATE spoa pthread sdsl suk slog ${ZLIB_LIBRARIES} ${HTS} OpenMP::OpenMP_CXX)
install(TARGETS hypo DESTINATION ${CMAKE_INSTALL_BINDIR})




