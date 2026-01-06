
#add_executable(forward forward.cpp)
#target_link_libraries(forward PRIVATE nlohmann_json::nlohmann_json)
#target_link_libraries(forward PUBLIC tomocam)
#target_link_libraries(forward PUBLIC TIFF::TIFF)

#add_executable(backproj backproj.cpp)
#target_link_libraries(backproj PRIVATE nlohmann_json::nlohmann_json)
#target_link_libraries(backproj PUBLIC tomocam)
#target_link_libraries(backproj PUBLIC TIFF::TIFF)

add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(test_fft PUBLIC tomocam)

add_executable(test_loss test_loss.cpp)
target_link_libraries(test_loss PUBLIC tomocam)

add_executable(test_gradient test_grad.cpp)
target_link_libraries(test_gradient PRIVATE tomocam)

add_executable(test_cg test_cg.cpp)
target_link_libraries(test_cg PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(test_cg PUBLIC tomocam)
target_link_libraries(test_cg PUBLIC TIFF::TIFF)
