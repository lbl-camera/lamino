
# Find FFTW library
#
# This module defines:
#   FFTW_FOUND
#   FFTW_INCLUDE_DIR
#   FFTW_LIBRARIES
#   FFTWF_LIBRARIES
#
# Targets:
#   FFTW::double   (double precision, links fftw3)
#   FFTW::float    (single precision, links fftw3f)
#
# Set FFTW_DIR to the directory containing FFTW3Config.cmake, e.g.:
#   cmake -DFFTW_DIR=~/fftw/lib/cmake/fftw3 ..
#
# For module-mode fallback set FFTW_ROOT or FFTW_DIR to the installation
# prefix (NERSC sets $FFTW_ROOT via environment modules).

include(FindPackageHandleStandardArgs)

# --- Config mode (preferred when available) ---
# Recent FFTW3 releases ship FFTW3Config.cmake; try the user-supplied dir first.
#
# Some vendor packages ship a broken config: the Cray cray-fftw module installs
# FFTW3Config.cmake but not the FFTW3LibraryDepends.cmake it include()s, and its
# paths still point at the RPM BUILDROOT. That include() is a hard error which
# CONFIG QUIET does not suppress, so probe the config for completeness before
# handing it to find_package() and fall through to module mode when it is broken.
function(_fftw_config_is_usable config_file out_var)
    set(${out_var} FALSE PARENT_SCOPE)
    if (NOT EXISTS "${config_file}")
        return()
    endif()
    get_filename_component(_dir "${config_file}" DIRECTORY)
    file(STRINGS "${config_file}" _includes REGEX "^[ \t]*include[ \t]*\\(")
    foreach (_line IN LISTS _includes)
        # Only relative-to-config includes are checkable without evaluating the file.
        if (_line MATCHES "\\$\\{CMAKE_CURRENT_LIST_DIR\\}/([^\"\\)]+)")
            if (NOT EXISTS "${_dir}/${CMAKE_MATCH_1}")
                message(STATUS
                    "Ignoring incomplete FFTW3 CMake package in ${_dir}: "
                    "missing ${CMAKE_MATCH_1}")
                return()
            endif()
        endif()
    endforeach()
    set(${out_var} TRUE PARENT_SCOPE)
endfunction()

# Locate a candidate config without loading it. FFTW_DIR may be either the
# config directory itself or an install prefix (NERSC's module sets it to
# <prefix>/lib), so search both directly and via the usual suffixes.
find_path(FFTW3_CONFIG_DIR
    NAMES
        FFTW3Config.cmake
        fftw3-config.cmake
    HINTS
        "${FFTW_DIR}"
        "$ENV{FFTW_DIR}"
        "${FFTW_ROOT}"
        "$ENV{FFTW_ROOT}"
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES
        cmake/fftw3
        lib/cmake/fftw3
        lib64/cmake/fftw3
        share/cmake/fftw3
)
mark_as_advanced(FFTW3_CONFIG_DIR)

set(_fftw3_usable FALSE)
foreach (_name FFTW3Config.cmake fftw3-config.cmake)
    if (FFTW3_CONFIG_DIR AND EXISTS "${FFTW3_CONFIG_DIR}/${_name}")
        _fftw_config_is_usable("${FFTW3_CONFIG_DIR}/${_name}" _fftw3_usable)
        break()
    endif()
endforeach()

if (_fftw3_usable)
    find_package(FFTW3 CONFIG QUIET
        PATHS
            "${FFTW3_CONFIG_DIR}"
        NO_DEFAULT_PATH
    )
endif()

if (FFTW3_FOUND)
    # Alias to our conventional target names if the config didn't create them.
    if (TARGET FFTW3::fftw3 AND NOT TARGET FFTW::double)
        add_library(FFTW::double ALIAS FFTW3::fftw3)
    endif()
    if (TARGET FFTW3::fftw3f AND NOT TARGET FFTW::float)
        add_library(FFTW::float ALIAS FFTW3::fftw3f)
    endif()
    message(STATUS "Found FFTW3 (config): ${FFTW3_DIR}")
    return()
endif()

# --- Module mode fallback ---
# FFTW_ROOT (CMake 3.12+) is honoured automatically by find_path/find_library.
# $ENV{FFTW_DIR} is the conventional HPC environment variable for the install prefix.
set(FFTW_SEARCH_PATHS
    /usr/local
    /opt/homebrew
    /opt/local
    /usr
    /opt
)

find_path(FFTW_INCLUDE_DIR
    NAMES
        fftw3.h
    HINTS
        $ENV{FFTW_ROOT}
        $ENV{FFTW_DIR}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES
        include
    PATHS
        ${FFTW_SEARCH_PATHS}
)

find_library(FFTW_LIBRARY
    NAMES
        fftw3
    HINTS
        $ENV{FFTW_ROOT}
        $ENV{FFTW_DIR}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES
        lib64
        lib
    PATHS
        ${FFTW_SEARCH_PATHS}
)

find_library(FFTWF_LIBRARY
    NAMES
        fftw3f
    HINTS
        $ENV{FFTW_ROOT}
        $ENV{FFTW_DIR}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES
        lib64
        lib
    PATHS
        ${FFTW_SEARCH_PATHS}
)

find_package_handle_standard_args(FFTW
    REQUIRED_VARS
        FFTW_INCLUDE_DIR
        FFTW_LIBRARY
        FFTWF_LIBRARY
)

if (FFTW_FOUND)
    set(FFTW_LIBRARIES  ${FFTW_LIBRARY})
    set(FFTWF_LIBRARIES ${FFTWF_LIBRARY})

    mark_as_advanced(
        FFTW_INCLUDE_DIR
        FFTW_LIBRARY
        FFTWF_LIBRARY
    )

    if (NOT TARGET FFTW::double)
        add_library(FFTW::double INTERFACE IMPORTED GLOBAL)
        target_include_directories(FFTW::double INTERFACE ${FFTW_INCLUDE_DIR})
        target_link_libraries(FFTW::double INTERFACE ${FFTW_LIBRARY})
    endif()
    message(STATUS "Found FFTW::double (module): ${FFTW_LIBRARY}")

    if (NOT TARGET FFTW::float)
        add_library(FFTW::float INTERFACE IMPORTED GLOBAL)
        target_include_directories(FFTW::float INTERFACE ${FFTW_INCLUDE_DIR})
        target_link_libraries(FFTW::float INTERFACE ${FFTWF_LIBRARY})
    endif()
    message(STATUS "Found FFTW::float (module): ${FFTWF_LIBRARY}")
else()
    message(FATAL_ERROR
        "FFTW not found.\n"
        "  Config mode:  set FFTW_DIR to the directory containing FFTW3Config.cmake\n"
        "                e.g. -DFFTW_DIR=~/fftw/lib/cmake/fftw3\n"
        "  Module mode:  set FFTW_ROOT to the installation prefix\n"
        "                e.g. -DFFTW_ROOT=~/fftw\n"
        "                or export FFTW_DIR=~/fftw (environment variable)"
    )
endif()
