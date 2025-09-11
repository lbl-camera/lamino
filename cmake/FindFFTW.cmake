
# Find or Fetch the FFTW library
#
# This module defines the following variables:
#   FFTW_FOUND
#   FFTW_INCLUDE_DIR
#   FFTW_LIBRARIES
#
# This module defines the following targets:
#   FFTW::fftw3
#   FFTW::fftw3f
# 



include(FindPackageHandleStandardArgs)

# try to find the local FFTW installation
set(FFTW_SEARCH_PATHS
    /usr/local
    /opt/local
    /opt/homebrew
    /usr
    /opt
)

find_path(FFTW_INCLUDE_DIR
    NAMES 
        fftw3.h
    HINTS 
        $ENV{FFTW_DIR} 
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        include
    PATHS 
        ${FFTW_SEARCH_PATHS}
)

find_library(FFTW_LIBRARIES
    NAMES 
        fftw3
    HINTS 
        $ENV{FFTW_DIR}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        lib64 
        lib
    PATHS 
        ${FFTW_SEARCH_PATHS}
)

find_library(FFTWF_LIBRARIES
    NAMES
        fftw3f
    HINTS
        $ENV{FFTW_DIR}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES 
        lib64 
        lib
    PATHS 
        ${FFTW_SEARCH_PATHS}
)

if (FFTW_INCLUDE_DIR AND FFTW_LIBRARIES AND FFTWF_LIBRARIES)
    set(FFTW_FOUND TRUE)


    FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW
        REQUIRED_VARS 
            FFTW_INCLUDE_DIR
            FFTW_LIBRARIES
            FFTWF_LIBRARIES
            FFTW_FOUND
    )

    mark_as_advanced(
        FFTW_INCLUDE_DIR
        FFTW_LIBRARIES
        FFTWF_LIBRARIES
    )

    # create imported target
    if (NOT TARGET FFTW::double)
        add_library(FFTW::double INTERFACE IMPORTED GLOBAL)
        target_include_directories(FFTW::double INTERFACE ${FFTW_INCLUDE_DIR})
        target_link_libraries(FFTW::double INTERFACE ${FFTW_LIBRARIES})
     endif()
     message(STATUS "Found FFTW::double ${FFTW_LIBRARIES}")

     if (NOT TARGET FFTW::float)
        add_library(FFTW::float INTERFACE IMPORTED GLOBAL)
        target_include_directories(FFTW::float INTERFACE ${FFTW_INCLUDE_DIR})
        target_link_libraries(FFTW::float INTERFACE ${FFTWF_LIBRARIES})
     endif()
     message(STATUS "Found FFTW::float ${FFTWF_LIBRARIES}")

else()
    message(FATAL_ERROR "local installation of FFTW not found")
endif()

