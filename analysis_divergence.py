"""Compute L2 norm of magnetic field divergence over time.

Expects out_bx_<step>.csv & out_by_<step>.csv inside Result/.
Outputs PNG.
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

l2=[]
for s in steps:
    bx = load(s,"bx").reshape(nx,ny).T
    by = load(s,"by").reshape(nx,ny).T
    div = (np.roll(bx,-1,1)-np.roll(bx,1,1))/(2*dx) + (np.roll(by,-1,0)-np.roll(by,1,0))/(2*dy)
    l2.append(np.sqrt(np.mean(div**2)))

plt.semilogy(steps,l2,'o-')
plt.xlabel("step"); plt.ylabel("L2(∇·B)")
plt.title("Magnetic Divergence Error")
plt.tight_layout()
plt.savefig("Result/divB_error.png",dpi=200)
print("Saved Result/divB_error.png")
