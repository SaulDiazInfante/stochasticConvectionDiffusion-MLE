import numpy as np
def read_fortran_binary_2d(filename ="../data/sampled_path.bin"):
    """
    Reads a 2D Fortran-written binary array (real64) using stream access.

    Parameters
    ----------
    filename : str
        Path to the binary file.
    shape : tuple of int
        Shape of the 2D array (rows, cols) as written in Fortran.

    Returns
    -------
    np.ndarray
        The 2D NumPy array in correct shape and order.
    """
    dtype = np.float64

    # Load data
    with open(filename, "rb") as f:
        shape = np.fromfile(f, dtype=np.int32, count=2)
        rows, cols = shape[0], shape[1]
        # Read the rest as float64 (Fortran real(real64))
        data = np.fromfile(f, dtype=np.float64)
        # Reshape with Fortran order (column-major)
        array = data.reshape((rows, cols), order='F')

    print("Array shape:", array.shape)
    return array
