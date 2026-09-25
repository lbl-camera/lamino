
add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft PUBLIC tomocam)
add_test(NAME test_fft COMMAND test_fft)

add_executable(test_cg test_cg.cpp)
target_link_libraries(test_cg PUBLIC tomocam)
target_link_libraries(test_cg PUBLIC TIFF::TIFF)
add_test(NAME test_cg COMMAND test_cg)

add_executable(test_fft_r2c test_fft_r2c.cpp)
target_link_libraries(test_fft_r2c PUBLIC tomocam)
add_test(NAME test_fft_r2c COMMAND test_fft_r2c)

add_executable(test_transpose test_transpose.cpp)
target_link_libraries(test_transpose PUBLIC tomocam)
add_test(NAME test_transpose COMMAND test_transpose)

add_executable(test_fftshift test_fftshift.cpp)
target_link_libraries(test_fftshift PUBLIC tomocam)
add_test(NAME test_fftshift COMMAND test_fftshift)

#add_executable(test_masking test_masking.cpp)
#target_link_libraries(test_masking PUBLIC tomocam)

add_executable(test_grad test_grad.cpp)
target_link_libraries(test_grad PUBLIC tomocam)
add_test(NAME test_grad COMMAND test_grad)

add_executable(test_adjoint test_adjoint.cpp)
target_link_libraries(test_adjoint PUBLIC tomocam)
add_test(NAME test_adjoint COMMAND test_adjoint)

add_executable(test_toeplitz test_toeplitz.cpp)
target_link_libraries(test_toeplitz PUBLIC tomocam)
add_test(NAME test_toeplitz COMMAND test_toeplitz)

if(CMAKE_CUDA_COMPILER)
    add_executable(test_gpu_conjgrad test_gpu_conjgrad.cu)
    target_link_libraries(test_gpu_conjgrad PUBLIC tomocam)
    add_test(NAME test_gpu_conjgrad COMMAND test_gpu_conjgrad)

    add_executable(test_gpu_gradient test_gpu_gradient.cu)
    target_link_libraries(test_gpu_gradient PUBLIC tomocam)
    add_test(NAME test_gpu_gradient COMMAND test_gpu_gradient)

    add_executable(test_gpu_alignment test_gpu_alignment.cu)
    target_link_libraries(test_gpu_alignment PUBLIC tomocam)
    add_test(NAME test_gpu_alignment COMMAND test_gpu_alignment)

    add_executable(test_gpu_toeplitz test_gpu_toeplitz.cu)
    target_link_libraries(test_gpu_toeplitz PUBLIC tomocam)
    add_test(NAME test_gpu_toeplitz COMMAND test_gpu_toeplitz)
endif()


add_executable(test_forward_consistency test_forward_consistency.cpp)
target_link_libraries(test_forward_consistency PUBLIC tomocam)
add_test(NAME test_forward_consistency COMMAND test_forward_consistency)
