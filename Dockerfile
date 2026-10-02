# syntax=docker/dockerfile:1
# GPU build:  docker build -t lamino .
#   --build-arg CUDA_ARCH=89   (GPU compute capability, no GPU needed at build time)
# CPU build:  docker build -t lamino-cpu --build-arg USE_CUDA=OFF .
# Run:        docker run --rm --gpus all -v $PWD:/data lamino recon <args>
ARG CUDA_VERSION=13.1.0
ARG UBUNTU=ubuntu24.04

FROM docker.io/nvidia/cuda:${CUDA_VERSION}-devel-${UBUNTU} AS build
ARG USE_CUDA=ON
ARG CUDA_ARCH=89
ARG MARCH=x86-64-v3
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake ninja-build git ca-certificates clang-18 libomp-18-dev \
        libtbb-dev libfftw3-dev libtiff-dev \
    && rm -rf /var/lib/apt/lists/*
ENV CC=clang-18 CXX=clang++-18

# FINUFFT (+ cuFINUFFT)
RUN git clone --depth 1 https://github.com/flatironinstitute/finufft.git /src/finufft \
    && cmake -S /src/finufft -B /src/finufft/build -G Ninja \
        -DCMAKE_INSTALL_PREFIX=/opt/finufft -DCMAKE_INSTALL_LIBDIR=lib64 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CUDA_HOST_COMPILER=clang++-18 \
        -DFINUFFT_USE_CUDA=${USE_CUDA} "-DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCH}" \
        -DFINUFFT_BUILD_TESTS=OFF -DBUILD_TESTING=OFF -DFINUFFT_BUILD_EXAMPLES=OFF \
        -DFINUFFT_BUILD_PYTHON=OFF -DFINUFFT_STATIC_LINKING=OFF \
    && cmake --build /src/finufft/build && cmake --install /src/finufft/build

# tomocam / recon / forward / adjoint
WORKDIR /src/lamino
COPY CMakeLists.txt ./
COPY cmake cmake
COPY include include
COPY src src
COPY tests/CMakeLists.txt tests/CMakeLists.txt
RUN cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/lamino \
        -DCMAKE_CUDA_HOST_COMPILER=clang++-18 \
        -DUSE_CUDA=${USE_CUDA} "-DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCH}" \
        -DTOMOCAM_MARCH=${MARCH} -DENABLE_TESTS=OFF \
        -Dfinufft_DIR=/opt/finufft/lib64/cmake/finufft \
    && cmake --build build && cmake --install build

FROM docker.io/nvidia/cuda:${CUDA_VERSION}-runtime-${UBUNTU}
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libomp5-18 libtbb12 libfftw3-single3 libfftw3-double3 libtiff6 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/finufft /opt/finufft
COPY --from=build /opt/lamino /opt/lamino
RUN printf '/opt/finufft/lib64\n/opt/lamino/lib\n' > /etc/ld.so.conf.d/lamino.conf && ldconfig
ENV PATH=/opt/lamino/bin:$PATH
WORKDIR /data
CMD ["recon"]
