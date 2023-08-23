import numpy as np
from memory_profiler import profile


@profile
def euclidean_broadcast(x, y):
    """Euclidean square distance matrix.

    Inputs:
    x: (N, m) numpy array
    y: (N, m) numpy array

    Ouput:
    (N, N) Euclidean square distance matrix:
    r_ij = (x_ij - y_ij)^2
    """
    diff = x[:, np.newaxis, :] - y[np.newaxis, :, :]

    return (diff * diff).sum(axis=2)


if __name__ == "__main__":
    nsamples = 2000
    nfeat = 50
    x = 10. * np.random.random([nsamples, nfeat])

    euclidean_broadcast(x, x)
