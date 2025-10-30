
find_package(nlohmann_json REQUIRED)
find_package(TIFF REQUIRED)

add_executable(forward forward.cpp)
target_link_libraries(forward PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(forward PUBLIC tomocam)
target_link_libraries(forward PUBLIC TIFF::TIFF)

add_executable(backproj backproj.cpp)
target_link_libraries(backproj PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(backproj PUBLIC tomocam)
target_link_libraries(backproj PUBLIC TIFF::TIFF)

add_executable(chessboard chessboard.cpp)
target_link_libraries(chessboard PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(chessboard PUBLIC tomocam)
target_link_libraries(chessboard PUBLIC TIFF::TIFF)

add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(test_fft PUBLIC tomocam)


add_executable(test_cg test_cg.cpp)
target_link_libraries(test_cg PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(test_cg PUBLIC tomocam)
target_link_libraries(test_cg PUBLIC TIFF::TIFF)

# add tests in order
add_test(NAME test01 COMMAND chessboard "${CMAKE_SOURCE_DIR}/config/chess01.json")
add_test(NAME test02 COMMAND chessboard "${CMAKE_SOURCE_DIR}/config/chess02.json")
add_test(NAME test03 COMMAND forward    "${CMAKE_SOURCE_DIR}/config/fwd01.json")
add_test(NAME test04 COMMAND forward    "${CMAKE_SOURCE_DIR}/config/fwd02.json")
add_test(NAME test05 COMMAND backproj   "${CMAKE_SOURCE_DIR}/config/bwd01.json")
