import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Grid parameters
nx, ny = 40, 40

# Create figure and axes
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Animation update function
def update(step):
    ax1.clear()
    ax2.clear()
    
    # Read pressure and velocity data
    data_p = pd.read_csv(f"output_p_{step}.csv", header=None)
    x = data_p[0].values.reshape(nx, ny)
    y = data_p[1].values.reshape(nx, ny)
    p = data_p[2].values.reshape(nx, ny)
    data_u = pd.read_csv(f"output_u_{step}.csv", header=None)
    data_v = pd.read_csv(f"output_v_{step}.csv", header=None)
    u = data_u[2].values.reshape(nx, ny)
    v = data_v[2].values.reshape(nx, ny)
    
    # Plot pressure field
    c1 = ax1.contourf(x, y, p, cmap="viridis", levels=20)
    ax1.set_title(f"Pressure Field (Step {step})")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect("equal")
    
    # Plot velocity field
    ax2.quiver(x[::2, ::2], y[::2, ::2], u[::2, ::2], v[::2, ::2])
    ax2.set_title(f"Velocity Field (Step {step})")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_aspect("equal")
    
    return c1, ax2

# Generate animation
steps = list(range(0, 50, 2))  # Steps: 0, 2, 4, ..., 48
ani = FuncAnimation(fig, update, frames=steps, blit=False)
plt.tight_layout()
ani.save("flow_animation.mp4", writer="ffmpeg", fps=10)
plt.show()