
add_executable(sirt main.cpp)
target_link_libraries(sirt PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(sirt PUBLIC tomocam)
target_link_libraries(sirt PUBLIC finufft::finufft)

add_executable(toeplitz toeplitz.cpp)
target_link_libraries(toeplitz PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(toeplitz PUBLIC tomocam)
target_link_libraries(toeplitz PUBLIC finufft::finufft)


add_executable(gradient test_grad.cpp)
target_link_libraries(gradient PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(gradient PUBLIC tomocam)
target_link_libraries(gradient PUBLIC finufft::finufft)


add_executable(ffttest test_fft.cpp)
target_link_libraries(ffttest PUBLIC tomocam)
target_link_libraries(ffttest PUBLIC fftw3)
