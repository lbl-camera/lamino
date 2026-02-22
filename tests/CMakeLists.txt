
# three-dimnsional shepp-log phantom
add_executable(shepp3d shepplogan.cpp)
target_link_libraries(shepp3d PUBLIC tomocam)
target_link_libraries(shepp3d PUBLIC TIFF::TIFF)
target_link_libraries(shepp3d PUBLIC tomlplusplus::tomlplusplus)

add_executable(forward forward.cpp)
target_link_libraries(forward PUBLIC tomocam)
target_link_libraries(forward PUBLIC TIFF::TIFF)
target_link_libraries(forward PUBLIC tomlplusplus::tomlplusplus)

add_executable(backproj backproj.cpp)
target_link_libraries(backproj PUBLIC tomocam)
target_link_libraries(backproj PUBLIC TIFF::TIFF)
target_link_libraries(backproj PUBLIC tomlplusplus::tomlplusplus)

add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft PUBLIC tomocam)

add_executable(test_gradient test_grad.cpp)
target_link_libraries(test_gradient PUBLIC tomocam)

add_executable(test_loss test_loss.cpp)
target_link_libraries(test_loss PUBLIC tomocam)

add_executable(test_lipschitz test_lipschitz.cpp)
target_link_libraries(test_lipschitz PUBLIC tomocam)


add_executable(test_cg test_cg.cpp)
target_link_libraries(test_cg PUBLIC tomocam)
target_link_libraries(test_cg PUBLIC TIFF::TIFF)
