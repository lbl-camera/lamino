
find_package(nlohmann_json REQUIRED)
find_package(TIFF REQUIRED)

add_executable(forward forward.cpp)
target_link_libraries(forward PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(forward PUBLIC tomocam)
target_link_libraries(forward PUBLIC TIFF::TIFF)
