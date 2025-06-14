import numpy as np

# Adjust this to match your Fortran array type
dtype = np.float64

# Load data
file_name ="../data/sampled_path.bin"

with open(file_name, "rb") as f:
    # Read two 4-byte integers (int32) for dimensions
    shape = np.fromfile(f, dtype=np.int32, count=2)
    rows, cols = shape[0], shape[1]

    # Read the rest as float64 (Fortran real(real64))
    data = np.fromfile(f, dtype=np.float64)

    # Reshape with Fortran order (column-major)
    array = data.reshape((rows, cols), order='F')

print("Array shape:", array.shape)


rows = data.size // cols
A = data.reshape((rows, cols))

print("Shape:", A.shape)
