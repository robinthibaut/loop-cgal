import loop_cgal
import numpy as np
def test():

    X,Y,Z = np.meshgrid(np.linspace(0, 1, 10),
                       np.linspace(0, 1, 10),
                          np.linspace(0, 1, 10))

    step_vector = np.array([.1, .1, .1])
    nsteps = np.ones(3)*10
    points, tri = loop_cgal.generate_mesh_from_numpy(Z,np.zeros(3),step_vector,nsteps,.6)
    print("tri", tri)
    print("points", points)
if __name__ == "__main__":
    test()
