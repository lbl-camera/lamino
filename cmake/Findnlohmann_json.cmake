# First try to find nlohmann_json via CONFIG mode (if installed with CMake)
find_package(nlohmann_json ${nlohmann_json_FIND_VERSION} CONFIG QUIET)

# If not found via CONFIG, fall back to MODULE mode
if(NOT nlohmann_json_FOUND)
    include(FindPackageHandleStandardArgs)

    # Search for the header file
    find_path(nlohmann_json_INCLUDE_DIR
        NAMES nlohmann/json.hpp
        HINTS 
            $ENV{nlohmann_json_DIR}
            ${nlohmann_json_ROOT}
        PATH_SUFFIXES 
            include
            include/nlohmann_json
        PATHS
            ~/json
            ~/nlohmann_json
            /usr/local
            /opt/local
            /usr
            /opt
    )

    # Handle standard find_package arguments
    find_package_handle_standard_args(nlohmann_json
        REQUIRED_VARS nlohmann_json_INCLUDE_DIR
        VERSION_VAR nlohmann_json_VERSION
    )

    # Mark cache variables as advanced
    mark_as_advanced(nlohmann_json_INCLUDE_DIR)

    # Create imported target if found and not already created
    if(nlohmann_json_FOUND AND NOT TARGET nlohmann_json::nlohmann_json)
        add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
        
        # Set include directories using the proper generator expression
        set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${nlohmann_json_INCLUDE_DIR}"
        )
        
        # Provide nlohmann_json_INCLUDE_DIRS for legacy compatibility
        set(nlohmann_json_INCLUDE_DIRS "${nlohmann_json_INCLUDE_DIR}")
    endif()
endif()
