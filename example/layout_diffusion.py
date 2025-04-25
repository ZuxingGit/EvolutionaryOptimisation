import numpy as np
import matplotlib.pyplot as plt

# Room initialization (room coordinates and function labels)
rooms = {
    "Kitchen": {"pos": np.array([0.1, 0.1]), "size": 0.2},
    "Living Room": {"pos": np.array([0.8, 0.8]), "size": 0.3},
    "Bedroom": {"pos": np.array([0.5, 0.5]), "size": 0.25},
    "Bathroom": {"pos": np.array([0.2, 0.7]), "size": 0.15},
}

# Target: functionally related rooms should be close
relations = {
    ("Kitchen", "Living Room"): 1.0,
    ("Living Room", "Bedroom"): 0.5,
    ("Bedroom", "Bathroom"): 0.3,
}


# Visualize initial layout
def plot_layout(rooms, title="Room Layout"):
    plt.figure(figsize=(8, 8))
    for name, data in rooms.items():
        plt.scatter(data["pos"][0], data["pos"][1], s=2000 * data["size"], label=name)
        plt.text(data["pos"][0], data["pos"][1], name, fontsize=12, ha="center")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True)
    plt.title(title)
    plt.legend()
    plt.show()


plot_layout(rooms, "Initial Layout")


# Diffusion process
def optimize_layout(rooms, relations, iterations=50, step_size=0.01):
    for _ in range(iterations):
        for (room1, room2), weight in relations.items():
            pos1, pos2 = rooms[room1]["pos"], rooms[room2]["pos"]
            direction = pos2 - pos1
            distance = np.linalg.norm(direction)
            if distance > 1e-3:  # in case of division by zero
                direction = direction / distance
                rooms[room1]["pos"] += step_size * weight * direction
                rooms[room2]["pos"] -= step_size * weight * direction
        # Clip positions to stay within the unit square
        for room in rooms.values():
            room["pos"] = np.clip(room["pos"], 0, 1)


# Optimize layout
optimize_layout(rooms, relations)
plot_layout(rooms, "Optimized Layout")
