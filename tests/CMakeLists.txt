
add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft PUBLIC tomocam)

add_executable(test_loss test_loss.cpp)
target_link_libraries(test_loss PUBLIC tomocam)

add_executable(test_gradient test_grad.cpp)
target_link_libraries(test_gradient PRIVATE tomocam)

add_executable(test_cg test_cg.cpp)
target_link_libraries(test_cg PUBLIC tomocam)
target_link_libraries(test_cg PUBLIC TIFF::TIFF)

add_executable(test_fft_r2c test_fft_r2c.cpp)
target_link_libraries(test_fft_r2c PUBLIC tomocam)

add_executable(test_transpose test_transpose.cpp)
target_link_libraries(test_transpose PUBLIC tomocam)

add_executable(test_demag test_demag.cpp)
target_link_libraries(test_demag PUBLIC tomocam)

