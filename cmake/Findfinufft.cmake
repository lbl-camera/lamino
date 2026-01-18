
# Find or Fetch the finufft library
#
# This module defines the following variables:
#   finufft_FOUND
#   finufft_INCLUDE_DIR
#   finufft_LIBRARIES
#
# This module defines the following targets:
#   finufft::finufft
#   finufft::cufinufft



include(FindPackageHandleStandardArgs)

# try to find the local finufft installation
set(finufft_SEARCH_PATHS
    $ENV{HOME}/finufft
    /usr/local
    /opt/homebrew
    /opt/local
    /usr
    /opt
)

find_path(finufft_INCLUDE_DIR
    NAMES 
        finufft.h 
    HINTS 
        $ENV{finufft_DIR} 
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        include
    PATHS 
        ${finufft_SEARCH_PATHS}
)

find_library(finufft_LIBRARY
    NAMES 
        finufft
    HINTS 
        $ENV{finufft_DIR} 
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        lib64 
        lib
    PATHS 
        ${finufft_SEARCH_PATHS}
)

find_library(cufinufft_LIBRARY
    NAMES 
        cufinufft
    HINTS 
        $ENV{finufft_DIR} 
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        lib64 
        lib
    PATHS 
        ${finufft_SEARCH_PATHS}
)

if (finufft_INCLUDE_DIR AND finufft_LIBRARY)
    set(finufft_FOUND TRUE)
    set(finufft_LIBRARIES ${finufft_LIBRARY})
    if (cufinufft_LIBRARY)
        list(APPEND finufft_LIBRARIES ${cufinufft_LIBRARY})
    endif()

    find_package_handle_standard_args(finufft 
        REQUIRED_VARS 
            finufft_INCLUDE_DIR
            finufft_LIBRARY
    )

    mark_as_advanced(
        finufft_INCLUDE_DIR
        finufft_LIBRARY
        cufinufft_LIBRARY
    )

    # create imported target for finufft
    if (NOT TARGET finufft::finufft)
        add_library(finufft::finufft UNKNOWN IMPORTED GLOBAL)
        set_target_properties(finufft::finufft PROPERTIES IMPORTED_LOCATION ${finufft_LIBRARY})
        target_include_directories(finufft::finufft INTERFACE ${finufft_INCLUDE_DIR})
    endif()

    # create imported target for cufinufft
    if (cufinufft_LIBRARY AND NOT TARGET finufft::cufinufft)
        add_library(finufft::cufinufft UNKNOWN IMPORTED GLOBAL)
        set_target_properties(finufft::cufinufft PROPERTIES IMPORTED_LOCATION ${cufinufft_LIBRARY})
        target_include_directories(finufft::cufinufft INTERFACE ${finufft_INCLUDE_DIR})
    endif()

    message(STATUS "Found finufft: ${finufft_LIBRARY}")
    if (cufinufft_LIBRARY)
        message(STATUS "Found cufinufft: ${cufinufft_LIBRARY}")
    endif()
else()
    message(FATAL_ERROR "local installation of finufft not found")
endif()
