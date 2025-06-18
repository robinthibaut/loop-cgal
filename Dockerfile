FROM quay.io/pypa/manylinux_2_28_x86_64

RUN dnf install -y \
    wget \
    cmake \
    make \
    gcc \
    gcc-c++ \
    gmp-devel \
    mpfr-devel \
    python3-pip

# Remove existing Boost package
RUN dnf remove -y boost-devel

# Install Boost 1.72.0 manually
RUN wget -q https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz && \
    tar -xzf boost_1_88_0.tar.gz && \
    cd boost_1_88_0 && \
    ./bootstrap.sh && \
    ./b2 install && \
    cd .. && rm -rf boost_1_88_0*

# Optional: Install Eigen manually
RUN wget -q https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && \
    tar -xzf eigen-3.4.0.tar.gz && \
    cd eigen-3.4.0 && mkdir build && cd build && \
    cmake .. && make install && \
    cd ../.. && rm -rf eigen-3.4.0*

# Optional: Install CGAL manually
RUN wget -q https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz && \
    tar -xf CGAL-6.0.1.tar.xz && \
    cd CGAL-6.0.1 && mkdir build && cd build && \
    cmake .. && make install && \
    cd ../.. && rm -rf CGAL-5.5.2*
