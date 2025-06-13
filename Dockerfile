FROM quay.io/pypa/manylinux_2_28_x86_64

RUN dnf install -y \
    wget \
    cmake \
    make \
    gcc \
    gcc-c++ \
    boost-devel \
    gmp-devel \
    mpfr-devel \
    python3-pip

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
