import numpy as np
import matplotlib.pyplot as plt

def plot_random_paths(A, num_paths=10, seed=None, title="Random Path Samples"):
    """
    Plot a random subset of paths from a 2D array A[i,j], where i indexes time and j indexes paths.

    Parameters
    ----------
    A : np.ndarray
        2D array with shape (n_obs, n_paths), where each column is a path
    num_paths : int
        Number of random paths to plot
    seed : int or None
        Random seed for reproducibility
    title : str
        Title of the plot for the random paths
    """
    n_obs, n_paths = A.shape

    if seed is not None:
        np.random.seed(seed)

    # Time axis
    delta = 0.00001
    t_horizon = delta * n_obs
    time = np.linspace(0, t_horizon, n_obs)

    # Random sample paths
    selected_indices = np.random.choice(n_paths, size=min(num_paths, n_paths), replace=False)
    plt.figure(figsize=(10, 6))
    for j in selected_indices:
        plt.plot(time, A[:, j], alpha=0.5)
    plt.xlabel("Time")
    plt.ylabel("Path value")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(loc='upper right', fontsize='small', ncol=2)

    # First sample path
    k = np.random.randint(0, 2499)
    plt.figure(figsize=(10, 4))
    plt.plot(time, A[:, k], color='black', label="Path 1")
    plt.xlabel("Time")
    plt.ylabel("Path value")
    plt.title(f"Sample Path ({k})")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()

    plt.show()
