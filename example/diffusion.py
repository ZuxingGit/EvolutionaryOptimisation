import numpy as np
import matplotlib.pyplot as plt

# Initialize layout
layout = np.random.rand(10, 2)  # Randomly generate positions for 10 rooms

# Target position (e.g., center point)
target = np.array([0.5, 0.5])

# Visualize initial layout
plt.scatter(layout[:, 0], layout[:, 1], c="blue", label="Initial")
plt.scatter(target[0], target[1], c="red", label="Target")

# Diffusion process
for i in range(10):  # Iterate 10 times
    layout += 0.1 * (target - layout)  # Move towards the target position

# Visualize final layout
plt.scatter(layout[:, 0], layout[:, 1], c="green", label="Final")
plt.legend()
plt.show()
