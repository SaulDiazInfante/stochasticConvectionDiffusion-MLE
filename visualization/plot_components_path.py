import numpy as np
from read_fortran_binary_2d import *
from plot_random_paths import *
import matplotlib.pyplot as plt

# Example usage
file_name = "../data/sampled_path.bin"
path = read_fortran_binary_2d(file_name)
print("path[0, 0] =", path[0, 0])
print("path.shape =", path.shape)

plot_random_paths(path, num_paths=1000, title="Random Path Samples")
plt.show()
