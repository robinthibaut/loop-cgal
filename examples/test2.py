import loop_cgal
import numpy as np


def test():

    X, Y, Z = np.meshgrid(
        np.linspace(0, 1, 10), np.linspace(0, 1, 10), np.linspace(0, 1, 10)
    )

    step_vector = np.array([0.1, 0.1, 0.1])
    nsteps = np.ones(3) * 10
    result = loop_cgal.compute_optimized_intersection(
        Z,X,0.5,.5,0,0,0,0.1,0.1,0.1
    )
    


if __name__ == "__main__":
    test()
