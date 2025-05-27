# Loop-cgal

Loop-cgal is a Python package for mesh processing operations using the  CGAL (Computational Geometry Algorithms Library). It is designed for efficient geometric computations using pyvista objects.

## Features

- Python bindings for CGAL using `pybind11`.
- Current features:
    - clipping of 3D triangular surfaces
- Future features:
    - Marching cubes algorithm for isosurface extraction
    - Boolean operations on marching cube meshes.

## Installation

### Prerequisites

- C++17 or later
- Python 3.11 or later
- CGAL library
- Boost library
- CMake 3.15 or later
- pybind11
- scikit-build
- pyvista

### Build and Install

1. Clone the repository:
   ```bash
   git clone https://github.com/Loop3D/loop-cgal.git
   cd loop-cgal
   pip install .
   ```
2. Alternatively, you can install it directly from PyPI:
   ```bash
   pip install loop-cgal
   ```