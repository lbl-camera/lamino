
add_executable(sirt main.cpp)
target_link_libraries(sirt PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(sirt PRIVATE tomocam)
target_link_libraries(sirt PRIVATE finufft::finufft)

add_executable(toeplitz toeplitz.cpp)
target_link_libraries(toeplitz PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(toeplitz PRIVATE tomocam)
target_link_libraries(toeplitz PRIVATE finufft::finufft)
