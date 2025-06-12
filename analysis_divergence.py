"""Analyse magnetic field divergence from simulation output.

This script reads ``out_bx_*.csv`` and ``out_by_*.csv`` files located
in ``Result/`` and computes both the L2 norm and the maximum absolute
value of ``∇·B`` for each available time step.  Two subplots showing the
evolution of these quantities are saved to ``Result/divB_error.png``.
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt, re, os, glob

files = glob.glob("Result/out_bx_*.csv")
if not files:
    raise RuntimeError("No B field output found. Did you recompile & rerun the solver?")

steps = sorted(int(re.findall(r"_bx_(\d+).csv",f)[0]) for f in files)

def load(step,prefix):
    d = pd.read_csv(f"Result/out_{prefix}_{step}.csv",header=None)
    return d.values[:,2]

def grid_shape():
    df = pd.read_csv(files[0],header=None)
    xs = np.unique(df[0]); ys = np.unique(df[1])
    return len(xs), len(ys), xs, ys

nx, ny, xs, ys = grid_shape()
dx = xs[1]-xs[0]; dy = ys[1]-ys[0]

l2 = []
max_abs = []
for s in steps:
    bx = load(s, "bx").reshape(nx, ny).T
    by = load(s, "by").reshape(nx, ny).T
    div = (
        (np.roll(bx, -1, 1) - np.roll(bx, 1, 1)) / (2 * dx)
        + (np.roll(by, -1, 0) - np.roll(by, 1, 0)) / (2 * dy)
    )
    l2.append(np.sqrt(np.mean(div ** 2)))
    max_abs.append(np.max(np.abs(div)))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.semilogy(steps, l2, 'o-')
ax1.set_xlabel('step')
ax1.set_ylabel('L2(∇·B)')
ax1.set_title('L2 Divergence')

ax2.semilogy(steps, max_abs, 'o-r')
ax2.set_xlabel('step')
ax2.set_ylabel('max |∇·B|')
ax2.set_title('Max Divergence')

fig.tight_layout()
fig.savefig('Result/divB_error.png', dpi=200)
print(f"Saved Result/divB_error.png. Final L2={l2[-1]:.3e}, max={max_abs[-1]:.3e}")

