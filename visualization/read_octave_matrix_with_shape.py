import numpy as np

def read_fortran_binary_matrix(filename, shape):
    """
    Reads a 2D matrix stored as raw float64 values in Fortran column-major order.

    Parameters
    ----------
    filename : str
        Path to the binary file.
    shape : tuple of int
        Shape of the 2D array (rows, cols)

    Returns
    -------
    np.ndarray
        Reconstructed NumPy array.
    """
    rows, cols = shape
    data = np.fromfile(filename, dtype=np.float64, count=rows * cols)
    return data.reshape((rows, cols), order='F')  # Fortran-style layout
